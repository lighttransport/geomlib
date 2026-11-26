/*
 * geomlib - Convex polygon intersection area
 * (see geomlib_polygon_intersection_area.h)
 *
 * Clean-room clip-and-shoelace. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_intersection_area.h"
#include "geomlib/geomlib_polygon_clip.h"

#include <math.h>
#include <string.h>

double geomlib_polygon_intersection_area(const geomlib_vec2 *a, size_t na,
                                         const geomlib_vec2 *b, size_t nb)
{
    geomlib_polygon2d clipped;
    double twice = 0.0;
    size_t i, n;

    if (!a || !b || na < 3 || nb < 3) return 0.0;
    memset(&clipped, 0, sizeof(clipped));
    if (geomlib_clip_polygon_convex(a, na, b, nb, &clipped) !=
        GEOMLIB_MESH_SUCCESS) return 0.0;
    n = clipped.vertex_count;
    for (i = 0; i < n; ++i) {
        geomlib_vec2 p = clipped.vertices[i];
        geomlib_vec2 q = clipped.vertices[(i + 1) % n];
        twice += (double)p.x * q.y - (double)q.x * p.y;
    }
    geomlib_polygon2d_free(&clipped);
    return 0.5 * fabs(twice);
}
