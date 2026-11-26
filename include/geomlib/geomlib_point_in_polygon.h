/*
 * geomlib - 2D point-in-polygon tests
 *
 * Ray-crossing containment and integer winding number for simple polygons.
 */

#ifndef GEOMLIB_POINT_IN_POLYGON_H
#define GEOMLIB_POINT_IN_POLYGON_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 1 when q is strictly inside, 0 outside. Boundary behaviour is undefined. */
int geomlib_point_in_polygon(const geomlib_vec2 *poly, size_t n, geomlib_vec2 q);

/* Integer winding number (nonzero means inside for simple polygons). */
int geomlib_polygon_winding(const geomlib_vec2 *poly, size_t n, geomlib_vec2 q);

/* Signed distance from q to the polygon boundary (positive outside when the
 * winding number is zero). */
float geomlib_polygon_signed_distance(const geomlib_vec2 *poly, size_t n,
                                      geomlib_vec2 q);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POINT_IN_POLYGON_H */
