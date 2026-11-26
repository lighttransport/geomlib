/*
 * geomlib - Capsule / triangle overlap (see geomlib_capsule_triangle.h)
 *
 * Clean-room closest-feature enumeration. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_capsule_triangle.h"
#include "geomlib/geomlib_montecarlo.h"
#include "geomlib/geomlib_segment_distance.h"
#include "geomlib/geomlib_segment_triangle.h"

#include <float.h>

static void consider(geomlib_segment_triangle_distance_result *best,
                     geomlib_vec3 pa, geomlib_vec3 pb)
{
    float d = geomlib_vec3_length(geomlib_vec3_sub(pa, pb));
    if (d < best->distance) {
        best->distance = d;
        best->point_a = pa;
        best->point_b = pb;
    }
}

int geomlib_segment_triangle_distance(
    geomlib_vec3 a, geomlib_vec3 b,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_segment_triangle_distance_result *out)
{
    geomlib_segment_triangle_distance_result r;
    geomlib_vec3 tri[3], edge_a, edge_b;
    int i;

    r.distance = FLT_MAX;
    r.point_a = geomlib_vec3_create(0, 0, 0);
    r.point_b = geomlib_vec3_create(0, 0, 0);
    if (!out) return 0;
    tri[0] = v0; tri[1] = v1; tri[2] = v2;
    {
        geomlib_segment_triangle_hit hit;
        if (geomlib_segment_triangle_intersect(a, b, v0, v1, v2, &hit)) {
            r.distance = 0.0f;
            r.point_a = hit.point;
            r.point_b = hit.point;
            *out = r;
            return 1;
        }
    }
    for (i = 0; i < 2; ++i) {
        geomlib_vec3 p = i == 0 ? a : b;
        consider(&r, p, geomlib_closest_point_triangle(p, v0, v1, v2));
    }
    for (i = 0; i < 3; ++i) {
        geomlib_segment_distance sd;
        edge_a = tri[i];
        edge_b = tri[(i + 1) % 3];
        if (geomlib_segment_closest_points(a, b, edge_a, edge_b, &sd) ==
            GEOMLIB_MESH_SUCCESS)
            consider(&r, sd.point_a, sd.point_b);
    }
    *out = r;
    return 1;
}

int geomlib_capsule_triangle_overlap(
    geomlib_vec3 a, geomlib_vec3 b, float r,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_capsule_triangle_result *out)
{
    geomlib_capsule_triangle_result res;
    geomlib_segment_triangle_distance_result d;

    res.intersects = 0;
    res.distance = 0.0f;
    res.point_a = geomlib_vec3_create(0, 0, 0);
    res.point_b = geomlib_vec3_create(0, 0, 0);
    if (!out || !(r >= 0.0f)) return 0;
    if (!geomlib_segment_triangle_distance(a, b, v0, v1, v2, &d)) return 0;
    res.point_a = d.point_a;
    res.point_b = d.point_b;
    res.distance = d.distance - r;
    res.intersects = (res.distance <= 0.0f);
    *out = res;
    return res.intersects;
}
