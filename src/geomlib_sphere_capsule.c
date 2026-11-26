/*
 * geomlib - Sphere / capsule overlap (see geomlib_sphere_capsule.h)
 *
 * Clean-room axis projection built on the capsule helper. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sphere_capsule.h"
#include "geomlib/geomlib_capsule.h"

int geomlib_sphere_capsule_overlap(geomlib_vec3 center, float radius,
                                   geomlib_vec3 a, geomlib_vec3 b,
                                   float cap_radius,
                                   geomlib_sphere_capsule_result *out)
{
    geomlib_sphere_capsule_result r;
    geomlib_vec3 axis_point;
    float d;

    r.intersects = 0;
    r.distance = 0.0f;
    r.point = geomlib_vec3_create(0, 0, 0);
    if (!out || !(radius >= 0.0f && cap_radius >= 0.0f)) return 0;
    d = geomlib_capsule_point_distance(center, a, b, 0.0f, &axis_point);
    r.point = axis_point;
    r.distance = d - (radius + cap_radius);
    r.intersects = (r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}
