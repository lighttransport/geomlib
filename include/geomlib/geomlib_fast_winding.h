/*
 * geomlib - Fast winding number queries
 *
 * Hierarchical far-field approximation of the generalized winding number of a
 * triangle mesh, accelerated by a bounding-volume tree of triangle clusters.
 */

#ifndef GEOMLIB_FAST_WINDING_H
#define GEOMLIB_FAST_WINDING_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct geomlib_winding_tree geomlib_winding_tree;

/* Build a winding-number acceleration tree. Returns NULL on failure. */
geomlib_winding_tree *geomlib_winding_tree_build(const geomlib_trimesh *mesh,
                                                 size_t leaf_size);
void geomlib_winding_tree_free(geomlib_winding_tree *tree);

/* Approximate winding number at a query point. */
float geomlib_fast_winding_number(const geomlib_winding_tree *tree,
                                  geomlib_vec3 point);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_FAST_WINDING_H */
