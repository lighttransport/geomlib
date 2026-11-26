/*
 * geomlib - Polyline arc-length tables and evaluation
 *
 * Builds a cumulative length table and evaluates points at a given arc length.
 */

#ifndef GEOMLIB_ARC_LENGTH_H
#define GEOMLIB_ARC_LENGTH_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Fill `cumulative` (n entries, cumulative[0] = 0) and optionally `total`.
 * Returns GEOMLIB_MESH_ERROR_INVALID for n < 1 or NULL input. */
int geomlib_polyline_arc_table(const geomlib_vec3 *points, size_t n,
                               double *cumulative, double *total);

/* Total length of the polyline. */
double geomlib_polyline_length(const geomlib_vec3 *points, size_t n);

/* Evaluate the polyline at arc length s, clamping to [0,total]. Requires the
 * table built by geomlib_polyline_arc_table. */
geomlib_vec3 geomlib_polyline_at_length(const geomlib_vec3 *points, size_t n,
                                        const double *cumulative, double s);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ARC_LENGTH_H */
