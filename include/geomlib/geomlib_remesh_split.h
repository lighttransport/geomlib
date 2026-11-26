/*
 * geomlib - Isotropic surface refinement by edge splitting
 *
 * Splits edges longer than a target into shorter triangles while preserving
 * the surface, boundary and manifoldness.
 */

#ifndef GEOMLIB_REMESH_SPLIT_H
#define GEOMLIB_REMESH_SPLIT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Refine `mesh` until its longest edge is <= max_edge_length (or the pass
 * budget runs out). Returns the number of splits performed. */
int geomlib_refine_edges(geomlib_trimesh *mesh, float max_edge_length,
                         int max_passes);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REMESH_SPLIT_H */
