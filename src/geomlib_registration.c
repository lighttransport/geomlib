/*
 * geomlib - Rigid ICP registration (see geomlib_registration.h)
 *
 * Iterated closest point with Kabsch rotation fitting (cross-covariance SVD
 * replaced by the closed-form quaternion eigenvalue method of Horn, which
 * only needs a symmetric 4x4 eigensolver already available in the sparse
 * module). C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_registration.h"
#include "geomlib/geomlib_sparse.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_icp_params geomlib_icp_params_default(void)
{
    geomlib_icp_params p;
    p.max_iterations = 50;
    p.tolerance = 1e-6f;
    p.max_correspondence_distance = 0.0f;
    p.estimate_scale = false;
    return p;
}

/* Best 3x3 rotation aligning centered point sets a -> b (Horn quaternion). */
static geomlib_mat3 fit_rotation(const geomlib_vec3 *a, const geomlib_vec3 *b, size_t n)
{
    double H[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double N[16] = {0};
    double eig[4], vec[16];
    geomlib_mat3 r = geomlib_mat3_identity();
    size_t i;
    double sxx, syy, szz, sxy, syx, sxz, szx, syz, szy;

    for (i = 0; i < n; ++i) {
        double ax=a[i].x, ay=a[i].y, az=a[i].z;
        double bx=b[i].x, by=b[i].y, bz=b[i].z;
        H[0][0]+=ax*bx; H[0][1]+=ax*by; H[0][2]+=ax*bz;
        H[1][0]+=ay*bx; H[1][1]+=ay*by; H[1][2]+=ay*bz;
        H[2][0]+=az*bx; H[2][1]+=az*by; H[2][2]+=az*bz;
    }
    sxx=H[0][0]; syy=H[1][1]; szz=H[2][2];
    sxy=H[0][1]; syx=H[1][0]; sxz=H[0][2]; szx=H[2][0];
    syz=H[1][2]; szy=H[2][1];

    /* Horn's symmetric 4x4 K matrix; eigenvector (x,y,z,w) maximizes the
     * rotation alignment. */
    N[0] = sxx+syy+szz;
    N[5] = sxx-syy-szz;
    N[10] = -sxx+syy-szz;
    N[15] = -sxx-syy+szz;
    N[1] = N[4] = syz-szy;
    N[2] = N[8] = szx-sxz;
    N[3] = N[12] = sxy-syx;
    N[6] = N[9] = sxy+syx;
    N[7] = N[13] = sxz+szx;
    N[11] = N[14] = syz+szy;

    if (geomlib_eigen_sym_dense(N, 4, eig, vec) != GEOMLIB_SPARSE_SUCCESS)
        return r;
    {
        /* Ascending eigenvalues: the largest is at index 3. */
        const double *v = vec + 12;
        double x=v[0], y=v[1], z=v[2], w=v[3];
        double norm = sqrt(w*w+x*x+y*y+z*z);
        if (norm < 1e-12) return r;
        x/=norm; y/=norm; z/=norm; w/=norm;
        /* Column-major: m[col*3 + row]. */
        r.m[0] = (float)(w*w+x*x-y*y-z*z);
        r.m[1] = (float)(2*(x*y+w*z));
        r.m[2] = (float)(2*(x*z-w*y));
        r.m[3] = (float)(2*(x*y-w*z));
        r.m[4] = (float)(w*w-x*x+y*y-z*z);
        r.m[5] = (float)(2*(y*z+w*x));
        r.m[6] = (float)(2*(x*z+w*y));
        r.m[7] = (float)(2*(y*z-w*x));
        r.m[8] = (float)(w*w-x*x-y*y+z*z);
    }
    return r;
}

float geomlib_icp_rigid(const geomlib_trimesh *source,
                        const geomlib_trimesh *target,
                        const geomlib_icp_params *params,
                        geomlib_mat4 *transform)
{
    geomlib_bvh bvh;
    geomlib_mat4 total = geomlib_mat4_identity();
    geomlib_vec3 *src = NULL;
    size_t n;
    int it, maxit;
    float tol, prev = INFINITY, rms = INFINITY;

    if (!source || !target || !transform || !source->vertices || !target->vertices)
        return INFINITY;
    n = source->vertex_count;
    if (n == 0 || target->triangle_count == 0) return INFINITY;
    maxit = params ? params->max_iterations : 50;
    tol = params ? params->tolerance : 1e-6f;

    src = (geomlib_vec3 *)malloc(n * sizeof(geomlib_vec3));
    if (!src) return INFINITY;
    memcpy(src, source->vertices, n * sizeof(geomlib_vec3));
    geomlib_bvh_build(&bvh, target, GEOMLIB_BVH_QUALITY_MEDIUM);

    for (it = 0; it < maxit; ++it) {
        geomlib_vec3 *a = NULL, *b = NULL;
        size_t count = 0, i;
        geomlib_mat3 r;
        geomlib_vec3 ca = {0,0,0}, cb = {0,0,0};
        a = (geomlib_vec3 *)malloc(n * sizeof(geomlib_vec3));
        b = (geomlib_vec3 *)malloc(n * sizeof(geomlib_vec3));
        if (!a || !b) { free(a); free(b); break; }
        for (i = 0; i < n; ++i) {
            geomlib_vec3 closest;
            float dist;
            if (geomlib_bvh_closest_point(&bvh, src[i], &closest, &dist, NULL)) {
                if (params && params->max_correspondence_distance > 0.0f &&
                    dist > params->max_correspondence_distance) continue;
                a[count] = src[i];
                b[count] = closest;
                ++count;
            }
        }
        if (count < 3) { free(a); free(b); break; }
        for (i = 0; i < count; ++i) {
            ca.x += a[i].x; ca.y += a[i].y; ca.z += a[i].z;
            cb.x += b[i].x; cb.y += b[i].y; cb.z += b[i].z;
        }
        ca = geomlib_vec3_mul(ca, 1.0f/(float)count);
        cb = geomlib_vec3_mul(cb, 1.0f/(float)count);
        for (i = 0; i < count; ++i) {
            a[i] = geomlib_vec3_sub(a[i], ca);
            b[i] = geomlib_vec3_sub(b[i], cb);
        }
        r = fit_rotation(a, b, count);
        free(a); free(b);

        /* Compose total = Translate(cb) * R * Translate(-ca) * total. */
        {
            geomlib_mat4 step = geomlib_mat4_identity();
            step.m[0]=r.m[0]; step.m[1]=r.m[1]; step.m[2]=r.m[2];
            step.m[4]=r.m[3]; step.m[5]=r.m[4]; step.m[6]=r.m[5];
            step.m[8]=r.m[6]; step.m[9]=r.m[7]; step.m[10]=r.m[8];
            step.m[12]=cb.x - (r.m[0]*ca.x + r.m[3]*ca.y + r.m[6]*ca.z);
            step.m[13]=cb.y - (r.m[1]*ca.x + r.m[4]*ca.y + r.m[7]*ca.z);
            step.m[14]=cb.z - (r.m[2]*ca.x + r.m[5]*ca.y + r.m[8]*ca.z);
            total = geomlib_mat4_mul(step, total);
        }
        /* Apply step to the working source. */
        {
            geomlib_mat4 step = geomlib_mat4_identity();
            step.m[0]=r.m[0]; step.m[1]=r.m[1]; step.m[2]=r.m[2];
            step.m[4]=r.m[3]; step.m[5]=r.m[4]; step.m[6]=r.m[5];
            step.m[8]=r.m[6]; step.m[9]=r.m[7]; step.m[10]=r.m[8];
            step.m[12]=cb.x - (r.m[0]*ca.x + r.m[3]*ca.y + r.m[6]*ca.z);
            step.m[13]=cb.y - (r.m[1]*ca.x + r.m[4]*ca.y + r.m[7]*ca.z);
            step.m[14]=cb.z - (r.m[2]*ca.x + r.m[5]*ca.y + r.m[8]*ca.z);
            for (size_t i = 0; i < n; ++i) {
                geomlib_vec4 p = geomlib_mat4_mul_vec4(step,
                    geomlib_vec4_create(src[i].x, src[i].y, src[i].z, 1.0f));
                src[i] = geomlib_vec3_create(p.x, p.y, p.z);
            }
        }
        {
            double e = 0.0;
            for (size_t i = 0; i < n; ++i) {
                geomlib_vec3 closest; float dist;
                if (geomlib_bvh_closest_point(&bvh, src[i], &closest, &dist, NULL))
                    e += dist*dist;
            }
            rms = (float)sqrt(e / n);
        }
        if (fabsf(prev - rms) < tol) break;
        prev = rms;
    }
    *transform = total;
    free(src);
    geomlib_bvh_free(&bvh);
    return rms;
}

void geomlib_transform_mesh(geomlib_trimesh *mesh, const geomlib_mat4 *m)
{
    size_t i;
    if (!mesh || !m || !mesh->vertices) return;
    for (i = 0; i < mesh->vertex_count; ++i) {
        geomlib_vec4 p = geomlib_mat4_mul_vec4(*m,
            geomlib_vec4_create(mesh->vertices[i].x, mesh->vertices[i].y,
                                mesh->vertices[i].z, 1.0f));
        mesh->vertices[i] = geomlib_vec3_create(p.x, p.y, p.z);
    }
}
