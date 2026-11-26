/*
 * geomlib - Sphere / AABB overlap (see geomlib_sphere_aabb.h)
 *
 * Clean-room clamp-and-compare distance test. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_sphere_aabb.h"

#include <math.h>

static float clampf(float v, float lo, float hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

int geomlib_sphere_aabb_overlap(geomlib_vec3 center, float radius,
                                geomlib_vec3 bmin, geomlib_vec3 bmax,
                                geomlib_sphere_aabb_result *out)
{
    geomlib_vec3 q;
    float dx, dy, dz, d;
    geomlib_sphere_aabb_result r;

    r.intersects = 0;
    r.closest = geomlib_vec3_create(0, 0, 0);
    r.distance = 0.0f;
    if (!out || !(radius >= 0.0f)) return 0;

    q.x = clampf(center.x, bmin.x, bmax.x);
    q.y = clampf(center.y, bmin.y, bmax.y);
    q.z = clampf(center.z, bmin.z, bmax.z);
    dx = center.x - q.x;
    dy = center.y - q.y;
    dz = center.z - q.z;
    d = sqrtf(dx * dx + dy * dy + dz * dz);

    r.closest = q;
    r.distance = d - radius;
    r.intersects = (r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}
