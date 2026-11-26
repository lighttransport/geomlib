/*
 * geomlib - Spatial sorting via Morton (Z-order) and Hilbert curves
 *
 * Maps integer grid coordinates to a 1D key preserving locality, so that
 * sorting keys clusters nearby points. Used to accelerate BVH construction
 * and cache-coherent traversal.
 */

#ifndef GEOMLIB_SPATIAL_SORT_H
#define GEOMLIB_SPATIAL_SORT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

uint64_t geomlib_morton2_encode(uint32_t x, uint32_t y);
uint64_t geomlib_morton3_encode(uint32_t x, uint32_t y, uint32_t z);
void geomlib_morton2_decode(uint64_t code, uint32_t *x, uint32_t *y);
void geomlib_morton3_decode(uint64_t code, uint32_t *x, uint32_t *y, uint32_t *z);

/* Hilbert order for a 2D grid of side `side` (a power of two, >= 2).
 * Coordinates must satisfy 0 <= x,y < side. */
uint64_t geomlib_hilbert2_encode(uint32_t side, uint32_t x, uint32_t y);
void geomlib_hilbert2_decode(uint32_t side, uint64_t d, uint32_t *x, uint32_t *y);

/* Return a permutation placing spatially close points first. `order` must
 * hold `count` entries. spatial_sort2 uses Hilbert, spatial_sort3 Morton. */
int geomlib_spatial_sort2(const geomlib_vec2 *points, size_t count,
                          uint32_t *order);
int geomlib_spatial_sort3(const geomlib_vec3 *points, size_t count,
                          uint32_t *order);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPATIAL_SORT_H */
