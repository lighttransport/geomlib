/*
 * geomlib - Polyline simplification (Douglas-Peucker)
 *
 * Keeps the endpoints and recursively inserts the point of maximum
 * perpendicular distance whenever it exceeds a tolerance.
 */

#ifndef GEOMLIB_SIMPLIFY_H
#define GEOMLIB_SIMPLIFY_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Simplify an open polyline. `indices` must hold `count` entries; on output
 * `out_count` gives the number of kept indices in original order. */
int geomlib_simplify_polyline(const geomlib_vec3 *points, size_t count,
                              float epsilon, uint32_t *indices,
                              size_t *out_count);

/* Maximum perpendicular distance from all input points to the kept polyline. */
float geomlib_simplify_max_error(const geomlib_vec3 *points, size_t count,
                                 const uint32_t *indices, size_t kept);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SIMPLIFY_H */
