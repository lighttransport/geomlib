/*
 * geomlib - 2D Poisson-disk sampling inside a polygon
 *
 * Bridson dart throwing constrained by a simple polygon domain.
 */

#ifndef GEOMLIB_POISSON2D_H
#define GEOMLIB_POISSON2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float radius;
    int candidates;
    uint32_t seed;
    size_t max_samples;
} geomlib_poisson2d_params;

/* Fill `out` (capacity `out_capacity`) with samples inside `polygon`. */
size_t geomlib_poisson_disk_polygon(const geomlib_vec2 *polygon, size_t np,
                                    const geomlib_poisson2d_params *params,
                                    geomlib_vec2 *out, size_t out_capacity);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POISSON2D_H */
