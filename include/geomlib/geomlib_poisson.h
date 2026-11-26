/*
 * geomlib - Poisson surface reconstruction
 *
 * Clean-room implementation of the Poisson method (Kazhdan et al. 2006):
 * oriented points define a vector field whose divergence is solved for an
 * indicator function, and the surface is extracted as an isosurface.
 */

#ifndef GEOMLIB_POISSON_H
#define GEOMLIB_POISSON_H

#include "geomlib_mesh.h"
#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t resolution;   /* Grid resolution on the longest axis (default 48) */
    float padding;         /* Extra padding as a fraction of bounds (default 0.1) */
    int max_iterations;    /* CG iterations (default 1000) */
    float tolerance;       /* CG relative tolerance (default 1e-6) */
} geomlib_poisson_params;

geomlib_poisson_params geomlib_poisson_params_default(void);

/*
 * Reconstruct a closed surface from an oriented point cloud. The cloud must
 * contain normals (geomlib_pointcloud_estimate_normals can generate them).
 * Returns GEOMLIB_MESH_SUCCESS and fills `out` on success.
 */
int geomlib_reconstruct_poisson(const geomlib_pointcloud *cloud,
                                const geomlib_poisson_params *params,
                                geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POISSON_H */
