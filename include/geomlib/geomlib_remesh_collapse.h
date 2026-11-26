/*
 * geomlib - Isotropic surface coarsening by edge collapse
 *
 * Collapses short interior edges when the link condition holds, preserving
 * manifoldness and optional boundaries.
 */

#ifndef GEOMLIB_REMESH_COLLAPSE_H
#define GEOMLIB_REMESH_COLLAPSE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Collapse edges shorter than `min_edge_length`. Returns the collapse count.
 * When `preserve_boundary` is nonzero, boundary edges are not collapsed. */
int geomlib_coarsen_edges(geomlib_trimesh *mesh, float min_edge_length,
                          int max_passes, int preserve_boundary);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REMESH_COLLAPSE_H */
