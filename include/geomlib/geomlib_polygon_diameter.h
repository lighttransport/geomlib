/*
 * geomlib - Diameter of a convex polygon
 *
 * Rotating-calipers O(n) computation of the farthest pair of vertices.
 */

#ifndef GEOMLIB_POLYGON_DIAMETER_H
#define GEOMLIB_POLYGON_DIAMETER_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Diameter of a convex polygon given in consistent (CCW or CW) order. Returns
 * 0 for fewer than two vertices. */
double geomlib_convex_polygon_diameter(const geomlib_vec2 *poly, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_DIAMETER_H */
