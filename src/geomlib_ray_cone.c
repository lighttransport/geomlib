/*
 * geomlib - Ray / cone intersection (see geomlib_ray_cone.h)
 *
 * Clean-room quadratic in cone coordinates. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_ray_cone.h"

#include <float.h>
#include <math.h>

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

int geomlib_ray_cone_intersect(geomlib_vec3 origin, geomlib_vec3 direction,
                               geomlib_vec3 apex, geomlib_vec3 axis,
                               float half_angle, geomlib_ray_cone_hit *out)
{
    double u[3], v[3], w[3], A, B, C, c, c2, wl, disc, sq, t0, t1, best;
    int i, found;

    if (!out || !(half_angle > 1e-4f && half_angle < 1.5707963f)) return 0;
    u[0] = (double)origin.x - apex.x;
    u[1] = (double)origin.y - apex.y;
    u[2] = (double)origin.z - apex.z;
    v[0] = direction.x; v[1] = direction.y; v[2] = direction.z;
    w[0] = axis.x; w[1] = axis.y; w[2] = axis.z;
    wl = sqrt(dot3(w, w));
    if (!(wl > 1e-12)) return 0;
    for (i = 0; i < 3; ++i) w[i] /= wl;
    c = cos(half_angle);
    c2 = c * c;
    A = pow(dot3(v, w), 2.0) - dot3(v, v) * c2;
    B = 2.0 * (dot3(u, w) * dot3(v, w) - dot3(u, v) * c2);
    C = pow(dot3(u, w), 2.0) - dot3(u, u) * c2;

    best = DBL_MAX;
    found = 0;
    if (fabs(A) < 1e-14) {
        if (fabs(B) > 1e-14) {
            double t = -C / B;
            double ax = dot3(u, w) + t * dot3(v, w);
            if (t >= 0.0 && ax >= 0.0 && t < best) { best = t; found = 1; }
        }
    } else {
        disc = B * B - 4.0 * A * C;
        if (disc >= 0.0) {
            sq = sqrt(disc);
            t0 = (-B - sq) / (2.0 * A);
            t1 = (-B + sq) / (2.0 * A);
            if (t0 > t1) { double t = t0; t0 = t1; t1 = t; }
            if (t0 >= 0.0 && dot3(u, w) + t0 * dot3(v, w) >= 0.0 &&
                t0 < best) { best = t0; found = 1; }
            if (t1 >= 0.0 && dot3(u, w) + t1 * dot3(v, w) >= 0.0 &&
                t1 < best) { best = t1; found = 1; }
        }
    }
    if (!found) return 0;
    out->intersects = 1;
    out->t = best;
    out->point = geomlib_vec3_create(
        (float)((double)direction.x * best + origin.x),
        (float)((double)direction.y * best + origin.y),
        (float)((double)direction.z * best + origin.z));
    return 1;
}
