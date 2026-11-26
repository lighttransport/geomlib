/*
 * geomlib - Anisotropic metric tensors (see geomlib_metric.h)
 *
 * Clean-room symmetric PSD eigen solve and curvature-driven construction.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_metric.h"
#include "geomlib/geomlib_curvature.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void outer(geomlib_vec3 a, geomlib_vec3 b, double o[3][3])
{
    double av[3] = {(double)a.x, (double)a.y, (double)a.z};
    double bv[3] = {(double)b.x, (double)b.y, (double)b.z};
    int i, j;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) o[i][j] = av[i] * bv[j];
}

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double o[3])
{
    o[0] = a[1] * b[2] - a[2] * b[1];
    o[1] = a[2] * b[0] - a[0] * b[2];
    o[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3(const double a[3])
{
    return sqrt(dot3(a, a));
}

geomlib_metric3 geomlib_metric_identity(void)
{
    geomlib_metric3 M;
    memset(M.m, 0, sizeof(M.m));
    M.m[0][0] = M.m[1][1] = M.m[2][2] = 1.0;
    return M;
}

geomlib_metric3 geomlib_metric_diagonal(double l0, double l1, double l2)
{
    geomlib_metric3 M = geomlib_metric_identity();
    M.m[0][0] = fabs(l0);
    M.m[1][1] = fabs(l1);
    M.m[2][2] = fabs(l2);
    return M;
}

geomlib_metric3 geomlib_metric_from_frame(const geomlib_vec3 axes[3],
                                          const double eig[3])
{
    geomlib_metric3 M;
    double o[3][3];
    int a, i, j;
    memset(M.m, 0, sizeof(M.m));
    for (a = 0; a < 3; ++a) {
        outer(axes[a], axes[a], o);
        for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                M.m[i][j] += fabs(eig[a]) * o[i][j];
    }
    for (i = 0; i < 3; ++i)
        for (j = i + 1; j < 3; ++j) M.m[j][i] = M.m[i][j];
    return M;
}

double geomlib_metric_norm(geomlib_metric3 M, geomlib_vec3 v)
{
    double x[3] = {(double)v.x, (double)v.y, (double)v.z};
    double s = 0.0;
    int i, j;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) s += x[i] * M.m[i][j] * x[j];
    return sqrt(fmax(0.0, s));
}

double geomlib_metric_distance(geomlib_metric3 M, geomlib_vec3 p,
                               geomlib_vec3 q)
{
    return geomlib_metric_norm(M, geomlib_vec3_sub(p, q));
}

double geomlib_metric_det(geomlib_metric3 M)
{
    const double *a = M.m[0], *b = M.m[1], *c = M.m[2];
    return a[0] * (b[1] * c[2] - b[2] * c[1])
         - a[1] * (b[0] * c[2] - b[2] * c[0])
         + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

/* Jacobi eigenvalues of a symmetric 3x3 matrix, descending. */
static void symmetric_eigen(double in[3][3], double eig[3])
{
    double a[3][3];
    int sweep, p, q;
    memcpy(a, in, sizeof(a));
    for (sweep = 0; sweep < 48; ++sweep) {
        double off = 0.0;
        for (p = 0; p < 2; ++p)
            for (q = p + 1; q < 3; ++q) off += fabs(a[p][q]);
        if (off < 1e-18) break;
        for (p = 0; p < 2; ++p)
            for (q = p + 1; q < 3; ++q) {
                double apq = a[p][q], tau, t, c, s;
                if (fabs(apq) < 1e-20) continue;
                tau = (a[q][q] - a[p][p]) / (2.0 * apq);
                t = (tau >= 0.0 ? 1.0 : -1.0) /
                    (fabs(tau) + sqrt(1.0 + tau * tau));
                c = 1.0 / sqrt(1.0 + t * t);
                s = t * c;
                {
                    double app = a[p][p], aqq = a[q][q];
                    int r;
                    a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                    a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                    a[p][q] = a[q][p] = 0.0;
                    for (r = 0; r < 3; ++r) {
                        if (r == p || r == q) continue;
                        {
                            double akp = a[r][p], akq = a[r][q];
                            a[r][p] = a[p][r] = c * akp - s * akq;
                            a[r][q] = a[q][r] = s * akp + c * akq;
                        }
                    }
                }
            }
    }
    eig[0] = a[0][0]; eig[1] = a[1][1]; eig[2] = a[2][2];
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (eig[j] > eig[i]) { double t = eig[i]; eig[i] = eig[j]; eig[j] = t; }
}

double geomlib_metric_anisotropy(geomlib_metric3 M)
{
    double eig[3];
    double a[3][3];
    memcpy(a, M.m, sizeof(a));
    symmetric_eigen(a, eig);
    if (!(eig[2] > 1e-18)) return 1.0;
    return sqrt(fmax(1.0, eig[0] / eig[2]));
}

geomlib_metric3 geomlib_metric_interpolate(geomlib_metric3 A,
                                           geomlib_metric3 B, double t)
{
    geomlib_metric3 R;
    int i, j;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            R.m[i][j] = (1.0 - t) * A.m[i][j] + t * B.m[i][j];
    return R;
}

static void orthonormal_frame(geomlib_vec3 normal, geomlib_vec3 d1,
                              geomlib_vec3 d2, double axes[3][3])
{
    double n[3] = {(double)normal.x, (double)normal.y, (double)normal.z};
    double u[3], v[3], len;
    len = norm3(n);
    if (len > 1e-12) { n[0] /= len; n[1] /= len; n[2] /= len; }
    else { n[0] = 0.0; n[1] = 0.0; n[2] = 1.0; }
    u[0] = (double)d1.x; u[1] = (double)d1.y; u[2] = (double)d1.z;
    if (norm3(u) < 1e-12) {
        u[0] = fabs(n[0]) < 0.9 ? 1.0 : 0.0;
        u[1] = fabs(n[0]) < 0.9 ? 0.0 : 1.0;
        u[2] = 0.0;
    }
    {
        double d = dot3(u, n);
        u[0] -= d * n[0]; u[1] -= d * n[1]; u[2] -= d * n[2];
    }
    len = norm3(u);
    if (len > 1e-12) { u[0] /= len; u[1] /= len; u[2] /= len; }
    else { u[0] = 1.0; u[1] = 0.0; u[2] = 0.0; }
    cross3(n, u, v);
    {
        double d = dot3(v, n);
        v[0] -= d * n[0]; v[1] -= d * n[1]; v[2] -= d * n[2];
    }
    len = norm3(v);
    if (len > 1e-12) { v[0] /= len; v[1] /= len; v[2] /= len; }
    axes[0][0] = u[0]; axes[0][1] = u[1]; axes[0][2] = u[2];
    axes[1][0] = v[0]; axes[1][1] = v[1]; axes[1][2] = v[2];
    axes[2][0] = n[0]; axes[2][1] = n[1]; axes[2][2] = n[2];
    (void)d2;
}

int geomlib_metric_from_curvature(geomlib_vec3 normal, float k1, float k2,
                                  geomlib_vec3 d1, geomlib_vec3 d2,
                                  double target_size, double normal_scale,
                                  geomlib_metric3 *out)
{
    double frame[3][3];
    geomlib_vec3 axisv[3];
    double eig[3], h1, h2, hn, kmax, kmin, hmin;
    int i;

    if (!out || !(target_size > 1e-9)) return 0;
    if (!(normal_scale > 1e-9)) normal_scale = 1.0;
    orthonormal_frame(normal, d1, d2, frame);
    kmax = fabs((double)k1);
    kmin = fabs((double)k2);
    if (kmax < kmin) { double t = kmax; kmax = kmin; kmin = t; }
    hmin = target_size * 0.1;
    h1 = target_size / sqrt(1.0 + kmax * target_size);
    h2 = target_size / sqrt(1.0 + kmin * target_size);
    hn = target_size * normal_scale;
    if (h1 < hmin) h1 = hmin;
    if (h2 < hmin) h2 = hmin;
    eig[0] = 1.0 / (h1 * h1);
    eig[1] = 1.0 / (h2 * h2);
    eig[2] = 1.0 / (hn * hn);
    for (i = 0; i < 3; ++i)
        axisv[i] = geomlib_vec3_create((float)frame[0][i],
                                      (float)frame[1][i],
                                      (float)frame[2][i]);
    *out = geomlib_metric_from_frame(axisv, eig);
    return 1;
}

int geomlib_metric_field_from_curvature(const geomlib_trimesh *mesh,
                                        double target_size,
                                        double normal_scale,
                                        geomlib_metric3 *out)
{
    geomlib_mesh_curvature curv;
    geomlib_vec3 *normals = NULL;
    size_t i;
    if (!mesh || !out || !mesh->vertices || !mesh->indices) return 0;
    {
        geomlib_trimesh tmp = *mesh;
        tmp.normals = NULL;
        geomlib_trimesh_compute_normals(&tmp);
        normals = tmp.normals;
    }
    curv = geomlib_curvature_create(mesh->vertex_count);
    if (!curv.curvatures) { free(normals); return 0; }
    if (geomlib_compute_curvature(mesh, &curv) != GEOMLIB_MESH_SUCCESS) {
        geomlib_curvature_free(&curv);
        free(normals);
        return 0;
    }
    for (i = 0; i < mesh->vertex_count; ++i) {
        geomlib_vertex_curvature *c = &curv.curvatures[i];
        geomlib_vec3 n = normals ? normals[i]
                                 : geomlib_vec3_create(0, 0, 1);
        if (!geomlib_metric_from_curvature(
                n, c->principal_max, c->principal_min,
                c->principal_dir_max, c->principal_dir_min,
                target_size, normal_scale, &out[i])) {
            geomlib_curvature_free(&curv);
            free(normals);
            return 0;
        }
    }
    geomlib_curvature_free(&curv);
    free(normals);
    return 1;
}
