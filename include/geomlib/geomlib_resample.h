/*
 * geomlib - Arc-length polyline resampling
 *
 * Redistributes points uniformly along a 3D polyline.
 */

#ifndef GEOMLIB_RESAMPLE_H
#define GEOMLIB_RESAMPLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Write `samples` equally spaced points along `points`, including both ends.
 * `out` must hold `samples` entries. */
int geomlib_resample_polyline(const geomlib_vec3 *points, size_t n,
                              size_t samples, geomlib_vec3 *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RESAMPLE_H */
