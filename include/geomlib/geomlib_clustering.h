/*
 * geomlib - Vertex-clustering mesh decimation
 *
 * Collapses vertices inside voxel cells to a single representative, then
 * rebuilds the triangle list, dropping degenerate and duplicate faces.
 */

#ifndef GEOMLIB_CLUSTERING_H
#define GEOMLIB_CLUSTERING_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Decimate a mesh by clustering vertices into voxels of `voxel_size`. */
int geomlib_decimate_vertex_clustering(const geomlib_trimesh *input,
                                       float voxel_size,
                                       geomlib_trimesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CLUSTERING_H */
