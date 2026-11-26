/*
 * geomlib - Dynamic BVH refit
 *
 * Updates an existing BVH's bounds after mesh vertices move, without
 * rebuilding the tree topology or primitive ordering.
 */

#ifndef GEOMLIB_BVH_REFIT_H
#define GEOMLIB_BVH_REFIT_H

#include "geomlib/geomlib_bvh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Refit BVH bounds to the current vertex positions of `mesh`. The mesh must
 * have the same topology as the one used to build the BVH. */
void geomlib_bvh_refit(geomlib_bvh *bvh, const geomlib_trimesh *mesh);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_BVH_REFIT_H */
