/*
 * geomlib - Minimum width of a convex polygon
 *
 * Rotating-calipers O(n) computation of the smallest distance between two
 * parallel supporting lines.
 */

#ifndef GEOMLIB_POLYGON_WIDTH_H
#define GEOMLIB_POLYGON_WIDTH_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Minimum width of a convex polygon in consistent (CCW or CW) order. Returns
 * 0 for fewer than three vertices. */
double geomlib_convex_polygon_min_width(const geomlib_vec2 *poly, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_WIDTH_H */
