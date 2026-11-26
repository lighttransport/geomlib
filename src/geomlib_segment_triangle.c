/*
 * geomlib - Segment / triangle intersection (see geomlib_segment_triangle.h)
 *
 * Clean-room Moller-Trumbore with a finite-parameter test. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_segment_triangle.h"

#include <float.h>
#include <math.h>
#include <string.h>

static void to_double(geomlib_vec3 v, double o[3])
{
    o[0] = (double)v.x; o[1] = (double)v.y; o[2] = (double)v.z;
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

int geomlib_segment_triangle_intersect(
    geomlib_vec3 a, geomlib_vec3 b,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_segment_triangle_hit *out)
{
    double A[3], B[3], V0[3], V1[3], V2[3];
    double dir[3], e1[3], e2[3], pvec[3], tvec[3], qvec[3];
    double det, invdet, u, v, t;
    const double eps = 1e-12;
    geomlib_segment_triangle_hit r;

    r.intersects = 0;
    r.t = 0.0f;
    r.point = geomlib_vec3_create(0, 0, 0);
    r.bary = geomlib_vec3_create(0, 0, 0);
    to_double(a, A); to_double(b, B);
    to_double(v0, V0); to_double(v1, V1); to_double(v2, V2);
    for (int i = 0; i < 3; ++i) dir[i] = B[i] - A[i];
    for (int i = 0; i < 3; ++i) e1[i] = V1[i] - V0[i];
    for (int i = 0; i < 3; ++i) e2[i] = V2[i] - V0[i];
    cross3(dir, e2, pvec);
    det = dot3(e1, pvec);
    if (fabs(det) < eps) return 0;   /* parallel or coplanar */
    invdet = 1.0 / det;
    for (int i = 0; i < 3; ++i) tvec[i] = A[i] - V0[i];
    u = dot3(tvec, pvec) * invdet;
    if (u < -1e-9 || u > 1.0 + 1e-9) return 0;
    cross3(tvec, e1, qvec);
    v = dot3(dir, qvec) * invdet;
    if (v < -1e-9 || u + v > 1.0 + 1e-9) return 0;
    t = dot3(e2, qvec) * invdet;
    if (t < -1e-9 || t > 1.0 + 1e-9) return 0;

    r.intersects = 1;
    r.t = (float)t;
    r.point = geomlib_vec3_create(
        (float)(A[0] + t * dir[0]),
        (float)(A[1] + t * dir[1]),
        (float)(A[2] + t * dir[2]));
    r.bary = geomlib_vec3_create((float)(1.0 - u - v), (float)u, (float)v);
    if (out) *out = r;
    return 1;
}
