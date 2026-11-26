/*
 * geomlib - Circle / AABB overlap (see geomlib_circle_aabb.h)
 *
 * Clean-room clamp-and-compare test. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_circle_aabb.h"

#include <math.h>

static float clampv(float v, float lo, float hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

int geomlib_circle_aabb_overlap(geomlib_vec2 center, float radius,
                                geomlib_vec2 bmin, geomlib_vec2 bmax,
                                geomlib_circle_aabb_result *out)
{
    geomlib_circle_aabb_result r;
    geomlib_vec2 q;
    float dx, dy, d;

    r.intersects = 0;
    r.distance = 0.0f;
    r.closest = geomlib_vec2_create(0, 0);
    if (!out || !(radius >= 0.0f)) return 0;
    q = geomlib_vec2_create(clampv(center.x, bmin.x, bmax.x),
                            clampv(center.y, bmin.y, bmax.y));
    dx = center.x - q.x;
    dy = center.y - q.y;
    d = sqrtf(dx * dx + dy * dy);
    r.closest = q;
    r.distance = d - radius;
    r.intersects = (r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}
