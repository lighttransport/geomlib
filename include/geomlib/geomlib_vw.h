/*
 * geomlib - Visvalingam-Whyatt polyline simplification
 *
 * Iteratively removes the interior point with the smallest effective area
 * (area of the triangle with its current neighbours), recomputing neighbours.
 */

#ifndef GEOMLIB_VW_H
#define GEOMLIB_VW_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Effective area of point `i` with its immediate neighbours (i must be an
 * interior point). Returns 0 for i<=0 or i>=count-1. */
double geomlib_vw_effective_area(const geomlib_vec3 *points, size_t count,
                                 size_t i);

/* Simplify an open polyline by repeatedly dropping the smallest-area interior
 * point while that area is <= area_epsilon. `indices` holds `count` entries;
 * `out_count` receives the number of kept indices in original order. */
int geomlib_simplify_polyline_vw(const geomlib_vec3 *points, size_t count,
                                 double area_epsilon, uint32_t *indices,
                                 size_t *out_count);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_VW_H */
