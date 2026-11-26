/*
 * geomlib - Sphere / triangle overlap (see geomlib_sphere_triangle.h)
 *
 * Clean-room closest-point test built on the existing triangle projection.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sphere_triangle.h"
#include "geomlib/geomlib_montecarlo.h"

#include <math.h>

int geomlib_sphere_triangle_overlap(
    geomlib_vec3 center, float radius,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_sphere_triangle_result *out)
{
    geomlib_sphere_triangle_result r;
    geomlib_vec3 closest;

    r.intersects = 0;
    r.closest = geomlib_vec3_create(0, 0, 0);
    r.distance = 0.0f;
    if (!out || !(radius >= 0.0f)) return 0;
    closest = geomlib_closest_point_triangle(center, v0, v1, v2);
    r.closest = closest;
    r.distance = geomlib_vec3_length(geomlib_vec3_sub(center, closest)) - radius;
    r.intersects = (r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}
