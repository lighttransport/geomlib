/*
 * geomlib - RANSAC shape detection (see geomlib_shape_detection.h)
 *
 * Clean-room implementation of the classic RANSAC paradigm (Fischler &
 * Bolles, 1981) for planes, spheres and cylinders. Hypothesis generation
 * from minimal samples, consensus counting by point-to-primitive distance,
 * and iterative refinement of the best model on its inliers.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_shape_detection.h"
#include "geomlib/geomlib_random.h"
#include "geomlib/geomlib_sparse.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_detected_shape geomlib_detected_shape_create(void)
{
    geomlib_detected_shape s;
    memset(&s, 0, sizeof(s));
    s.threshold = 0.01f;
    return s;
}

void geomlib_detected_shape_free(geomlib_detected_shape *s)
{
    if (!s) return;
    free(s->inliers);
    s->inliers = NULL;
    s->inlier_count = 0;
}

static float dist_to_plane(geomlib_vec3 p, geomlib_vec3 q, geomlib_vec3 n)
{
    return fabsf(geomlib_vec3_dot(geomlib_vec3_sub(p, q), n));
}

static float dist_to_sphere(geomlib_vec3 p, geomlib_vec3 c, float r)
{
    return fabsf(geomlib_vec3_length(geomlib_vec3_sub(p, c)) - r);
}

/* Distance to an infinite cylinder through c with unit axis d. */
static float dist_to_cylinder(geomlib_vec3 p, geomlib_vec3 c, geomlib_vec3 d, float r)
{
    geomlib_vec3 v = geomlib_vec3_sub(p, c);
    float t = geomlib_vec3_dot(v, d);
    geomlib_vec3 proj = geomlib_vec3_sub(v, geomlib_vec3_mul(d, t));
    return fabsf(geomlib_vec3_length(proj) - r);
}

/* Refine a plane's point and normal with PCA-like eigenvector iteration. */
static void refine_plane(const geomlib_pointcloud *cloud, const uint32_t *idx,
                         size_t n, geomlib_vec3 *point, geomlib_vec3 *normal)
{
    double cov[9] = {0};
    geomlib_vec3 c = {0,0,0};
    size_t i, j;
    for (i = 0; i < n; ++i) {
        geomlib_vec3 p = cloud->points[idx[i]];
        c.x += p.x; c.y += p.y; c.z += p.z;
    }
    c = geomlib_vec3_mul(c, 1.0f/(float)n);
    for (i = 0; i < n; ++i) {
        geomlib_vec3 d = geomlib_vec3_sub(cloud->points[idx[i]], c);
        double a[3] = {d.x, d.y, d.z};
        for (j = 0; j < 3; ++j)
            for (size_t k = 0; k < 3; ++k) cov[j*3+k] += a[j]*a[k];
    }
    /* Smallest eigenvector of the covariance is the plane normal; use the
     * symmetric eigensolver instead of an ad-hoc power iteration. */
    {
        double dense[9] = {
            cov[0], cov[1], cov[2],
            cov[3], cov[4], cov[5],
            cov[6], cov[7], cov[8]
        };
        double eig[3], vec[9];
        if (geomlib_eigen_sym_dense(dense, 3, eig, vec) != GEOMLIB_SPARSE_SUCCESS)
            return;
        *point = c;
        *normal = geomlib_vec3_create((float)vec[0], (float)vec[1], (float)vec[2]);
    }
}

/* Algebraic sphere fit from inliers (linear least squares on x^2+y^2+z^2). */
static void refine_sphere(const geomlib_pointcloud *cloud, const uint32_t *idx,
                          size_t n, geomlib_vec3 *center, float *radius)
{
    double A[4][4] = {{0}}, b[4] = {0};
    size_t i, j, k;
    for (i = 0; i < n; ++i) {
        geomlib_vec3 p = cloud->points[idx[i]];
        double x=p.x, y=p.y, z=p.z;
        double row[4] = {x,y,z,1.0};
        double rhs = x*x+y*y+z*z;
        for (j = 0; j < 4; ++j) {
            b[j] += row[j]*rhs;
            for (k = 0; k < 4; ++k) A[j][k] += row[j]*row[k];
        }
    }
    /* Gauss elimination with partial pivoting. */
    for (k = 0; k < 4; ++k) {
        size_t pivot = k;
        for (i = k+1; i < 4; ++i) if (fabs(A[i][k]) > fabs(A[pivot][k])) pivot = i;
        if (fabs(A[pivot][k]) < 1e-12) return;
        if (pivot != k) {
            for (j = 0; j < 4; ++j) { double t=A[k][j]; A[k][j]=A[pivot][j]; A[pivot][j]=t; }
            double t=b[k]; b[k]=b[pivot]; b[pivot]=t;
        }
        for (i = k+1; i < 4; ++i) {
            double f = A[i][k]/A[k][k];
            for (j = k; j < 4; ++j) A[i][j] -= f*A[k][j];
            b[i] -= f*b[k];
        }
    }
    {
        double x[4];
        for (k = 4; k-- > 0;) {
            double s = b[k];
            for (j = k+1; j < 4; ++j) s -= A[k][j]*x[j];
            x[k] = s/A[k][k];
        }
        *center = geomlib_vec3_create((float)(x[0]*0.5), (float)(x[1]*0.5), (float)(x[2]*0.5));
        *radius = (float)sqrt(x[3] + x[0]*x[0]*0.25 + x[1]*x[1]*0.25 + x[2]*x[2]*0.25);
    }
}

geomlib_detected_shape *geomlib_detect_shapes_ransac(
    const geomlib_pointcloud *cloud,
    geomlib_shape_type type,
    float threshold,
    size_t max_shapes,
    uint32_t seed,
    size_t *out_count)
{
    geomlib_pcg32 rng;
    bool *available;
    geomlib_detected_shape *result;
    size_t found = 0;
    size_t i, iter, n;
    const size_t iterations = 300;

    if (!cloud || !out_count || !cloud->points || cloud->point_count == 0)
        return NULL;
    if (!(threshold > 0.0f)) return NULL;
    n = cloud->point_count;
    available = (bool *)malloc(n * sizeof(bool));
    result = (geomlib_detected_shape *)calloc(max_shapes ? max_shapes : 1,
                                              sizeof(geomlib_detected_shape));
    if (!available || !result) {
        free(available); free(result);
        return NULL;
    }
    for (i = 0; i < n; ++i) available[i] = true;
    geomlib_pcg32_init(&rng, seed ? seed : 0x1234abcdu, 0x9e3779b9u);

    while (found < max_shapes) {
        size_t best_count = 0;
        uint32_t *best_inliers = NULL;
        geomlib_vec3 best_center = {0,0,0}, best_normal = {0,0,1};
        float best_radius = 0.0f;
        size_t remaining = 0;
        for (i = 0; i < n; ++i) if (available[i]) remaining++;
        if (remaining < 3) break;

        for (iter = 0; iter < iterations; ++iter) {
            uint32_t pick[3], count = 0, attempts = 0;
            while (count < 3 && attempts++ < 100) {
                uint32_t p = (uint32_t)(geomlib_pcg32_nextf(&rng) * (float)n);
                if (available[p]) {
                    bool dup = false;
                    for (i = 0; i < count; ++i) if (pick[i] == p) dup = true;
                    if (!dup) pick[count++] = p;
                }
            }
            if (count < 3) continue;

            geomlib_vec3 center = {0,0,0}, normal = {0,0,1};
            float radius = 0.0f;
            if (type == GEOMLIB_SHAPE_PLANE) {
                geomlib_vec3 a = cloud->points[pick[0]], b = cloud->points[pick[1]],
                             c = cloud->points[pick[2]];
                geomlib_vec3 nrm = geomlib_vec3_cross(geomlib_vec3_sub(b,a),
                                                      geomlib_vec3_sub(c,a));
                if (geomlib_vec3_length(nrm) < 1e-8f) continue;
                center = a;
                normal = geomlib_vec3_normalize(nrm);
            } else if (type == GEOMLIB_SHAPE_SPHERE) {
                /* Circumcenter of three points in 3D is underdetermined; use
                 * midpoint of the farthest pair as a seed and let the
                 * algebraic fit refine it. */
                geomlib_vec3 a = cloud->points[pick[0]], b = cloud->points[pick[1]];
                center = geomlib_vec3_mul(geomlib_vec3_add(a,b), 0.5f);
                radius = geomlib_vec3_length(geomlib_vec3_sub(a,b)) * 0.5f;
                if (!(radius > 0.0f)) continue;
            } else {
                /* Cylinder: axis from a sampled point direction, refined later. */
                geomlib_vec3 a = cloud->points[pick[0]], b = cloud->points[pick[1]];
                geomlib_vec3 axis = geomlib_vec3_sub(b, a);
                if (geomlib_vec3_length(axis) < 1e-8f) continue;
                center = a;
                normal = geomlib_vec3_normalize(axis);
                /* Seed radius from the third point's distance to the axis. */
                {
                    geomlib_vec3 v = geomlib_vec3_sub(cloud->points[pick[2]], center);
                    float t = geomlib_vec3_dot(v, normal);
                    radius = geomlib_vec3_length(
                        geomlib_vec3_sub(v, geomlib_vec3_mul(normal, t)));
                }
                if (!(radius > 0.0f)) continue;
            }

            {
                uint32_t *inliers = (uint32_t *)malloc(n * sizeof(uint32_t));
                size_t cnt = 0;
                if (!inliers) { free(best_inliers); free(inliers); goto fail; }
                for (i = 0; i < n; ++i) {
                    if (!available[i]) continue;
                    float d;
                    if (type == GEOMLIB_SHAPE_PLANE)
                        d = dist_to_plane(cloud->points[i], center, normal);
                    else if (type == GEOMLIB_SHAPE_SPHERE)
                        d = dist_to_sphere(cloud->points[i], center, radius);
                    else
                        d = dist_to_cylinder(cloud->points[i], center, normal, radius);
                    if (d <= threshold) inliers[cnt++] = (uint32_t)i;
                }
                if (cnt > best_count) {
                    free(best_inliers);
                    best_inliers = inliers;
                    best_count = cnt;
                    best_center = center;
                    best_normal = normal;
                    best_radius = radius;
                } else {
                    free(inliers);
                }
            }
        }

        if (best_count == 0 || best_count < 3) { free(best_inliers); break; }

        /* Refine the winning model on its consensus set. */
        if (type == GEOMLIB_SHAPE_PLANE)
            refine_plane(cloud, best_inliers, best_count, &best_center, &best_normal);
        else if (type == GEOMLIB_SHAPE_SPHERE)
            refine_sphere(cloud, best_inliers, best_count, &best_center, &best_radius);
        else {
            /* Keep the sampled cylinder model; refinement is future work. */
        }

        {
            geomlib_detected_shape *s = &result[found];
            *s = geomlib_detected_shape_create();
            s->type = type;
            s->threshold = threshold;
            s->center = best_center;
            s->normal = best_normal;
            s->radius = best_radius;
            s->inlier_count = best_count;
            s->inliers = best_inliers;
            best_inliers = NULL;
            for (i = 0; i < best_count; ++i) available[s->inliers[i]] = false;
        }
        found++;
    }

    free(available);
    *out_count = found;
    return result;

fail:
    free(available);
    for (i = 0; i < found; ++i) geomlib_detected_shape_free(&result[i]);
    free(result);
    return NULL;
}

void geomlib_detected_shapes_free(geomlib_detected_shape *shapes, size_t count)
{
    size_t i;
    if (!shapes) return;
    for (i = 0; i < count; ++i) geomlib_detected_shape_free(&shapes[i]);
    free(shapes);
}
