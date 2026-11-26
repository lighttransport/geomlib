/*
 * geomlib - K-d tree for 3D point sets
 *
 * Balanced k-d tree with k-nearest-neighbour and radius queries.
 */

#ifndef GEOMLIB_KDTREE_H
#define GEOMLIB_KDTREE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct geomlib_kdtree geomlib_kdtree;

geomlib_kdtree *geomlib_kdtree_build(const geomlib_vec3 *points, size_t count);
void geomlib_kdtree_free(geomlib_kdtree *tree);

/* Find up to k nearest point indices, sorted by distance. Returns count. */
size_t geomlib_kdtree_knn(const geomlib_kdtree *tree, geomlib_vec3 query,
                          size_t k, uint32_t *indices, float *distances);

/* Find all points within radius. Returns count (up to max_results). */
size_t geomlib_kdtree_radius(const geomlib_kdtree *tree, geomlib_vec3 query,
                             float radius, uint32_t *indices,
                             float *distances, size_t max_results);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_KDTREE_H */
