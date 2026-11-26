/*
 * geomlib - RBF variational surface reconstruction
 *
 * Fits a compactly supported radial basis function implicit field to oriented
 * points and extracts its zero level set with marching cubes.
 */

#ifndef GEOMLIB_RBF_RECONSTRUCT_H
#define GEOMLIB_RBF_RECONSTRUCT_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t resolution;   /* Grid resolution on the longest axis */
    float support;          /* RBF support radius (0 = auto from spacing) */
    int max_points;         /* Cap on sampled input points (default 400) */
} geomlib_rbf_params;

geomlib_rbf_params geomlib_rbf_params_default(void);

int geomlib_reconstruct_rbf(const geomlib_pointcloud *cloud,
                            const geomlib_rbf_params *params,
                            geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RBF_RECONSTRUCT_H */
