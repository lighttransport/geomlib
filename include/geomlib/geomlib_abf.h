/*
 * geomlib - Angle-based flattening (ABF-style)
 *
 * Minimizes the squared difference between 3D and 2D triangle corner angles
 * by gradient descent, starting from a free-boundary LSCM parameterization.
 */

#ifndef GEOMLIB_ABF_H
#define GEOMLIB_ABF_H

#include "geomlib_parameterization.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int iterations;      /* Gradient descent steps (default 60) */
    float step;          /* Learning rate (default 0.05) */
    float tolerance;     /* Early stop on max angle error (default 1e-3) */
} geomlib_abf_params;

geomlib_abf_params geomlib_abf_params_default(void);

/* Flatten `mesh` while pinning vertices `pin0` and `pin1`. */
int geomlib_parameterize_abf(const geomlib_trimesh *mesh,
                             uint32_t pin0, uint32_t pin1,
                             const geomlib_abf_params *params,
                             geomlib_uv_coords *out);

/* Maximum absolute corner-angle error (radians) between the 3D mesh and 2D UVs. */
float geomlib_abf_angle_error(const geomlib_trimesh *mesh,
                              const geomlib_uv_coords *uvs);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ABF_H */
