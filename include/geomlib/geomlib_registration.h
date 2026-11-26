/*
 * geomlib - Rigid registration (ICP)
 *
 * Clean-room implementation of the iterated closest point algorithm with
 * Kabsch/horn quaternion rotation fitting on paired points.
 */

#ifndef GEOMLIB_REGISTRATION_H
#define GEOMLIB_REGISTRATION_H

#include "geomlib_mesh.h"
#include "geomlib_bvh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int max_iterations;     /* default 50 */
    float tolerance;        /* RMS improvement tolerance, default 1e-6 */
    float max_correspondence_distance; /* 0 = unlimited */
    bool estimate_scale;    /* uniform scale (default false) */
} geomlib_icp_params;

geomlib_icp_params geomlib_icp_params_default(void);

/*
 * Align `source` onto `target`, writing the rigid transform as a 4x4
 * column-major matrix (applied to source points). Returns the final RMS error.
 */
float geomlib_icp_rigid(const geomlib_trimesh *source,
                        const geomlib_trimesh *target,
                        const geomlib_icp_params *params,
                        geomlib_mat4 *transform);

/* Transform a mesh in place by a column-major 4x4 matrix. */
void geomlib_transform_mesh(geomlib_trimesh *mesh, const geomlib_mat4 *m);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REGISTRATION_H */
