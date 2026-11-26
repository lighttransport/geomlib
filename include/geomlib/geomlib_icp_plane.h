/*
 * geomlib - Point-to-plane ICP registration (Chen & Medioni)
 *
 * Fits a rigid transform by minimizing normal distances between transformed
 * source vertices and their closest points on the target surface.
 */

#ifndef GEOMLIB_ICP_PLANE_H
#define GEOMLIB_ICP_PLANE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int max_iterations;   /* default 30 */
    float tolerance;      /* RMS convergence tolerance (default 1e-4) */
} geomlib_icp_plane_params;

geomlib_icp_plane_params geomlib_icp_plane_params_default(void);

/* Fit a transform placing `source` onto `target`. Writes the column-major
 * rigid transform and the final normal-distance RMS. */
int geomlib_icp_point_to_plane(const geomlib_trimesh *source,
                               const geomlib_trimesh *target,
                               const geomlib_icp_plane_params *params,
                               geomlib_mat4 *transform, float *final_rms);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ICP_PLANE_H */
