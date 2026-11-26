/*
 * geomlib - Triangle / AABB overlap (see geomlib_tri_aabb.h)
 *
 * Clean-room separating-axis test using the 13 candidate axes. C11, no
 * external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_tri_aabb.h"

#include <float.h>
#include <math.h>

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

static double norm3(const double a[3])
{
    return sqrt(dot3(a, a));
}

/* Project the (box-centred) triangle and the box half-extents onto an axis. */
static int separated(double v[3][3], double h[3], double axis[3])
{
    double mn = DBL_MAX, mx = -DBL_MAX, r, p;
    int i;
    for (i = 0; i < 3; ++i) {
        p = dot3(v[i], axis);
        if (p < mn) mn = p;
        if (p > mx) mx = p;
    }
    r = h[0] * fabs(axis[0]) + h[1] * fabs(axis[1]) + h[2] * fabs(axis[2]);
    return mn > r + 1e-12 || mx < -r - 1e-12;
}

int geomlib_triangle_aabb_overlap(geomlib_vec3 v0, geomlib_vec3 v1,
                                  geomlib_vec3 v2,
                                  geomlib_vec3 bmin, geomlib_vec3 bmax)
{
    double tri[3][3], h[3], center[3], e[3][3], n[3], axis[3];
    double bmin_a[3], bmax_a[3];
    double basis[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    int i, j;

    to_double(v0, tri[0]); to_double(v1, tri[1]); to_double(v2, tri[2]);
    to_double(bmin, bmin_a);
    to_double(bmax, bmax_a);
    for (i = 0; i < 3; ++i) center[i] = 0.5 * (bmin_a[i] + bmax_a[i]);
    for (i = 0; i < 3; ++i) h[i] = 0.5 * fabs(bmax_a[i] - bmin_a[i]);
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) tri[i][j] -= center[j];
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) e[i][j] = tri[(i + 1) % 3][j] - tri[i][j];

    /* Box face normals. */
    for (i = 0; i < 3; ++i)
        if (separated(tri, h, basis[i])) return 0;

    /* Triangle normal. */
    cross3(e[0], e[1], n);
    if (norm3(n) > 1e-18 && separated(tri, h, n)) return 0;

    /* Cross products of box axes with triangle edges. */
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            cross3(basis[i], e[j], axis);
            if (norm3(axis) > 1e-18 && separated(tri, h, axis)) return 0;
        }
    return 1;
}
