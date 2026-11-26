/*
 * geomlib - Point-cloud octree
 *
 * Adaptive octree over a point set with leaf capacity splitting, supporting
 * radius and k-nearest-neighbour queries.
 */

#ifndef GEOMLIB_OCTREE_H
#define GEOMLIB_OCTREE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct geomlib_octree geomlib_octree;

typedef struct {
    size_t leaf_capacity;   /* Points per leaf before splitting (default 8) */
    int max_depth;          /* Maximum subdivision depth (default 12) */
} geomlib_octree_params;

geomlib_octree_params geomlib_octree_params_default(void);

geomlib_octree *geomlib_octree_build(const geomlib_vec3 *points, size_t count,
                                     const geomlib_octree_params *params);
void geomlib_octree_free(geomlib_octree *tree);

size_t geomlib_octree_radius(const geomlib_octree *tree, geomlib_vec3 query,
                             float radius, uint32_t *indices,
                             float *distances, size_t max_results);
size_t geomlib_octree_knn(const geomlib_octree *tree, geomlib_vec3 query,
                          size_t k, uint32_t *indices, float *distances);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_OCTREE_H */
