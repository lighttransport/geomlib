/*
 * geomlib - Plane / plane intersection (see geomlib_plane_plane.h)
 *
 * Clean-room cross-product solution of the two plane equations. C11, no
 * external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_plane_plane.h"

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

int geomlib_plane_plane_intersection(geomlib_vec3 p0, geomlib_vec3 n0,
                                     geomlib_vec3 p1, geomlib_vec3 n1,
                                     geomlib_plane_line *out)
{
    double a[3], b[3], dir[3], m[3], q[3];
    double d0, d1, len, len2;
    int i;

    if (!out) return 0;
    memset(out, 0, sizeof(*out));
    to_double(n0, a);
    to_double(n1, b);
    len = sqrt(dot3(a, a));
    if (!(len > 1e-12)) return 0;
    for (i = 0; i < 3; ++i) a[i] /= len;
    len = sqrt(dot3(b, b));
    if (!(len > 1e-12)) return 0;
    for (i = 0; i < 3; ++i) b[i] /= len;

    cross3(a, b, dir);
    len2 = dot3(dir, dir);
    if (len2 < 1e-16) { out->parallel = 1; return 0; }
    len = sqrt(len2);

    {
        double pp0[3], pp1[3];
        to_double(p0, pp0);
        to_double(p1, pp1);
        d0 = dot3(a, pp0);
        d1 = dot3(b, pp1);
    }
    for (i = 0; i < 3; ++i) m[i] = d0 * b[i] - d1 * a[i];
    cross3(m, dir, q);
    for (i = 0; i < 3; ++i) q[i] /= len2;
    for (i = 0; i < 3; ++i) dir[i] /= len;

    out->point = geomlib_vec3_create((float)q[0], (float)q[1], (float)q[2]);
    out->direction = geomlib_vec3_create((float)dir[0], (float)dir[1],
                                         (float)dir[2]);
    return 1;
}
