/*
 * geomlib - Intersection area of two convex polygons
 *
 * Clips the first polygon against the second and applies the shoelace formula.
 */

#ifndef GEOMLIB_POLYGON_INTERSECTION_AREA_H
#define GEOMLIB_POLYGON_INTERSECTION_AREA_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

double geomlib_polygon_intersection_area(const geomlib_vec2 *a, size_t na,
                                         const geomlib_vec2 *b, size_t nb);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_INTERSECTION_AREA_H */
