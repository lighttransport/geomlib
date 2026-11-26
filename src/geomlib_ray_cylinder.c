/*
 * geomlib - Ray / cylinder intersection (see geomlib_ray_cylinder.h)
 *
 * Clean-room quadratic in the plane perpendicular to the axis. C11, no
 * external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_ray_cylinder.h"

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

int geomlib_ray_cylinder_intersect(
    geomlib_vec3 origin, geomlib_vec3 direction,
    geomlib_vec3 axis_point, geomlib_vec3 axis_dir, float radius,
    geomlib_ray_cylinder_hit *out)
{
    double O[3], D[3], P[3], W[3], v[3], dp[3];
    double a, b, c, disc, sq, t, t0, t1, axial, len;
    int i;
    geomlib_ray_cylinder_hit r;

    r.t = 0.0;
    r.point = geomlib_vec3_create(0, 0, 0);
    r.normal = geomlib_vec3_create(0, 0, 0);
    if (!out || !(radius >= 0.0f)) return 0;
    to_double(origin, O); to_double(direction, D);
    to_double(axis_point, P); to_double(axis_dir, W);
    len = sqrt(dot3(W, W));
    if (!(len > 1e-12)) return 0;
    for (i = 0; i < 3; ++i) W[i] /= len;
    for (i = 0; i < 3; ++i) v[i] = O[i] - P[i];
    axial = dot3(v, W);
    for (i = 0; i < 3; ++i) v[i] -= axial * W[i];
    axial = dot3(D, W);
    for (i = 0; i < 3; ++i) dp[i] = D[i] - axial * W[i];

    a = dot3(dp, dp);
    if (a < 1e-18) return 0;   /* parallel to the axis */
    b = dot3(v, dp);
    c = dot3(v, v) - (double)radius * radius;
    disc = b * b - a * c;
    if (disc < 0.0) return 0;
    sq = sqrt(disc);
    t0 = (-b - sq) / a;
    t1 = (-b + sq) / a;
    if (t1 < 0.0) return 0;
    t = (t0 >= 0.0) ? t0 : t1;

    r.t = t;
    r.point = geomlib_vec3_create(
        (float)(O[0] + t * D[0]),
        (float)(O[1] + t * D[1]),
        (float)(O[2] + t * D[2]));
    {
        double hit[3], n[3];
        hit[0] = (double)r.point.x; hit[1] = (double)r.point.y;
        hit[2] = (double)r.point.z;
        for (i = 0; i < 3; ++i) n[i] = hit[i] - P[i];
        axial = dot3(n, W);
        for (i = 0; i < 3; ++i) n[i] -= axial * W[i];
        len = sqrt(dot3(n, n));
        if (len > 1e-12) {
            for (i = 0; i < 3; ++i) n[i] /= len;
            r.normal = geomlib_vec3_create((float)n[0], (float)n[1], (float)n[2]);
        }
    }
    *out = r;
    return 1;
}
