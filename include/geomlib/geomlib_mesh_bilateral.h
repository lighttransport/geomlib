/*
 * geomlib - Feature-preserving bilateral mesh smoothing
 *
 * Averages neighbour positions with weights that decay with spatial distance
 * and normal disagreement, so sharp creases are preserved.
 */

#ifndef GEOMLIB_MESH_BILATERAL_H
#define GEOMLIB_MESH_BILATERAL_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float sigma_spatial;   /* Spatial Gaussian width (0 = mean edge length) */
    float sigma_normal;    /* Normal-similarity width (default 0.35) */
    int neighbors;         /* k nearest neighbours (default 12) */
    int iterations;        /* Smoothing passes (default 2) */
} geomlib_mesh_bilateral_params;

geomlib_mesh_bilateral_params geomlib_mesh_bilateral_params_default(void);

/* Smooth mesh vertex positions in place. Vertex normals are computed internally. */
int geomlib_mesh_bilateral_smooth(geomlib_trimesh *mesh,
                                  const geomlib_mesh_bilateral_params *params);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MESH_BILATERAL_H */
