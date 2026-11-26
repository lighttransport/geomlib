/*
 * geomlib - Moving least squares (MLS) surface projection
 *
 * Projects a query point onto the local weighted planar approximation of a
 * point cloud, using k-nearest neighbours and a symmetric 3x3 eigen solve.
 */

#ifndef GEOMLIB_MLS_H
#define GEOMLIB_MLS_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int neighbors;      /* k nearest neighbours (default 16) */
    int iterations;     /* Projection iterations (default 3) */
    float sigma;        /* Gaussian width (0 = local point spacing) */
} geomlib_mls_params;

geomlib_mls_params geomlib_mls_params_default(void);

/* Project `query` onto the MLS surface. On success writes the projected point
 * and, when non-NULL, the local normal. Returns GEOMLIB_MESH_SUCCESS. */
int geomlib_mls_project(const geomlib_pointcloud *cloud, geomlib_vec3 query,
                        const geomlib_mls_params *params,
                        geomlib_vec3 *projected, geomlib_vec3 *normal);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MLS_H */
