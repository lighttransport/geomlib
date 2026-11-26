/*
 * geomlib - Capsule collision (see geomlib_capsule.h)
 *
 * Built on the existing finite segment closest-point solver. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_capsule.h"
#include "geomlib/geomlib_segment_distance.h"

#include <math.h>

static float clamp01(float v)
{
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}

int geomlib_capsule_capsule_overlap(geomlib_vec3 a0, geomlib_vec3 a1, float r0,
                                    geomlib_vec3 b0, geomlib_vec3 b1, float r1,
                                    geomlib_capsule_result *out)
{
    geomlib_capsule_result r;
    geomlib_segment_distance sd;

    r.intersects = 0;
    r.distance = 0.0f;
    r.point_a = geomlib_vec3_create(0, 0, 0);
    r.point_b = geomlib_vec3_create(0, 0, 0);
    if (!out || !(r0 >= 0.0f && r1 >= 0.0f)) return 0;
    if (geomlib_segment_closest_points(a0, a1, b0, b1, &sd) != GEOMLIB_MESH_SUCCESS)
        return 0;
    r.point_a = sd.point_a;
    r.point_b = sd.point_b;
    r.distance = sd.distance - (r0 + r1);
    r.intersects = (r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}

float geomlib_capsule_point_distance(geomlib_vec3 q, geomlib_vec3 a0,
                                     geomlib_vec3 a1, float r,
                                     geomlib_vec3 *axis_point)
{
    geomlib_vec3 u = geomlib_vec3_sub(a1, a0);
    float u2 = geomlib_vec3_dot(u, u);
    float t = 0.0f;
    geomlib_vec3 closest;
    if (u2 > 1e-12f)
        t = clamp01(geomlib_vec3_dot(geomlib_vec3_sub(q, a0), u) / u2);
    closest = geomlib_vec3_add(a0, geomlib_vec3_mul(u, t));
    if (axis_point) *axis_point = closest;
    return geomlib_vec3_length(geomlib_vec3_sub(q, closest)) - r;
}
