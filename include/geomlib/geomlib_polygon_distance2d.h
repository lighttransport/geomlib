/*
 * geomlib - Distance between two simple 2D polygons
 *
 * Returns zero for overlapping polygons, otherwise the minimum distance
 * between their boundaries with witness points.
 */

#ifndef GEOMLIB_POLYGON_DISTANCE2D_H
#define GEOMLIB_POLYGON_DISTANCE2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float distance;
    geomlib_vec2 point_a;
    geomlib_vec2 point_b;
} geomlib_polygon_distance2d_result;

int geomlib_polygon_distance2d(const geomlib_vec2 *a, size_t na,
                               const geomlib_vec2 *b, size_t nb,
                               geomlib_polygon_distance2d_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_DISTANCE2D_H */
