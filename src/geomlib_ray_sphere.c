/*
 * geomlib - Ray / sphere intersection (see geomlib_ray_sphere.h)
 *
 * Clean-room analytic quadratic. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_ray_sphere.h"

#include <math.h>

int geomlib_ray_sphere_intersect(geomlib_vec3 origin, geomlib_vec3 direction,
                                 geomlib_vec3 center, float radius,
                                 geomlib_ray_sphere_hit *out)
{
    double ox = (double)origin.x - center.x;
    double oy = (double)origin.y - center.y;
    double oz = (double)origin.z - center.z;
    double dx = direction.x, dy = direction.y, dz = direction.z;
    double a = dx * dx + dy * dy + dz * dz;
    double b, c, disc, sq, t0, t1;

    if (!out || !(radius >= 0.0f) || !(a > 1e-18)) return 0;
    b = 2.0 * (dx * ox + dy * oy + dz * oz);
    c = ox * ox + oy * oy + oz * oz - (double)radius * radius;
    disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return 0;
    sq = sqrt(disc);
    t0 = (-b - sq) / (2.0 * a);
    t1 = (-b + sq) / (2.0 * a);
    if (t0 > t1) { double t = t0; t0 = t1; t1 = t; }

    out->t0 = t0;
    out->t1 = t1;
    out->p0 = geomlib_vec3_create((float)(origin.x + t0 * dx),
                                  (float)(origin.y + t0 * dy),
                                  (float)(origin.z + t0 * dz));
    out->p1 = geomlib_vec3_create((float)(origin.x + t1 * dx),
                                  (float)(origin.y + t1 * dy),
                                  (float)(origin.z + t1 * dz));
    {
        geomlib_vec3 n0 = geomlib_vec3_sub(out->p0, center);
        geomlib_vec3 n1 = geomlib_vec3_sub(out->p1, center);
        float l0 = geomlib_vec3_length(n0);
        float l1 = geomlib_vec3_length(n1);
        out->n0 = l0 > 1e-12f ? geomlib_vec3_mul(n0, 1.0f / l0)
                              : geomlib_vec3_create(0, 0, 0);
        out->n1 = l1 > 1e-12f ? geomlib_vec3_mul(n1, 1.0f / l1)
                              : geomlib_vec3_create(0, 0, 0);
    }
    return 1;
}
