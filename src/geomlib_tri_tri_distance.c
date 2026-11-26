/*
 * geomlib - Triangle / triangle distance (see geomlib_tri_tri_distance.h)
 *
 * Clean-room enumeration of vertex-triangle and edge-edge candidates. C11,
 * no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_tri_tri_distance.h"
#include "geomlib/geomlib_tri_tri.h"
#include "geomlib/geomlib_montecarlo.h"
#include "geomlib/geomlib_segment_distance.h"

#include <float.h>

static void consider(geomlib_tri_tri_distance_result *best,
                     geomlib_vec3 pa, geomlib_vec3 pb)
{
    float d = geomlib_vec3_length(geomlib_vec3_sub(pa, pb));
    if (d < best->distance) {
        best->distance = d;
        best->point_a = pa;
        best->point_b = pb;
    }
}

int geomlib_triangle_triangle_distance(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2,
    geomlib_vec3 b0, geomlib_vec3 b1, geomlib_vec3 b2,
    geomlib_tri_tri_distance_result *out)
{
    geomlib_tri_tri_distance_result r;
    geomlib_vec3 A[3], B[3];
    int i, j;

    r.distance = FLT_MAX;
    r.point_a = geomlib_vec3_create(0, 0, 0);
    r.point_b = geomlib_vec3_create(0, 0, 0);
    if (!out) return 0;
    A[0] = a0; A[1] = a1; A[2] = a2;
    B[0] = b0; B[1] = b1; B[2] = b2;

    if (geomlib_triangle_intersection(a0, a1, a2, b0, b1, b2).intersects) {
        r.distance = 0.0f;
        *out = r;
        return 1;
    }

    for (i = 0; i < 3; ++i) {
        consider(&r, A[i], geomlib_closest_point_triangle(A[i], B[0], B[1], B[2]));
        consider(&r, geomlib_closest_point_triangle(B[i], A[0], A[1], A[2]), B[i]);
    }
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            geomlib_segment_distance sd;
            if (geomlib_segment_closest_points(A[i], A[(i + 1) % 3],
                                               B[j], B[(j + 1) % 3],
                                               &sd) == GEOMLIB_MESH_SUCCESS)
                consider(&r, sd.point_a, sd.point_b);
        }
    *out = r;
    return 1;
}
