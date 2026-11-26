/*
 * geomlib - Iterative Laplacian curve smoothing
 *
 * Smooths an open polyline by averaged-neighbour displacements while
 * preserving the endpoints.
 */

#ifndef GEOMLIB_CURVE_SMOOTH_H
#define GEOMLIB_CURVE_SMOOTH_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Smooth `points` in place for `iterations` passes with weight `lambda`. */
int geomlib_smooth_polyline(geomlib_vec3 *points, size_t n,
                             int iterations, float lambda);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CURVE_SMOOTH_H */
