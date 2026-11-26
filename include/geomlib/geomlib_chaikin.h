/*
 * geomlib - Chaikin corner-cutting polyline subdivision
 *
 * Repeatedly replaces each segment with two points at 1/4 and 3/4 along it,
 * smoothing the curve while preserving its endpoints.
 */

#ifndef GEOMLIB_CHAIKIN_H
#define GEOMLIB_CHAIKIN_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Subdivide a polyline `iterations` times. `out` holds up to
 * `(count-1)*2^iterations + 1` points; returns the written count or 0. */
size_t geomlib_chaikin_subdivide(const geomlib_vec3 *points, size_t count,
                                 int iterations, geomlib_vec3 *out,
                                 size_t out_capacity);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CHAIKIN_H */
