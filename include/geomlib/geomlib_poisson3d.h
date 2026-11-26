/*
 * geomlib - 3D volume Poisson-disk sampling
 *
 * Bridson's dart-throwing algorithm with a background grid for acceleration,
 * generating a blue-noise point set inside an axis-aligned box.
 */

#ifndef GEOMLIB_POISSON3D_H
#define GEOMLIB_POISSON3D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 box_min, box_max;
    float radius;        /* minimum separation */
    int candidates;      /* attempts per active sample (default 30) */
    uint32_t seed;       /* RNG seed (0 = fixed default) */
    size_t max_samples;
} geomlib_poisson3d_params;

/* Fill `out` (capacity `out_capacity`) with samples; returns count written. */
size_t geomlib_poisson_disk_volume(const geomlib_poisson3d_params *params,
                                   geomlib_vec3 *out, size_t out_capacity);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POISSON3D_H */
