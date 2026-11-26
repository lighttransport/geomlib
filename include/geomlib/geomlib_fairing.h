/*
 * geomlib - Implicit mean-curvature flow / fairing
 *
 * Solves (M + lambda * L) x_new = M x_old with cotangent stiffness L and
 * lumped mass M, an unconditionally stable mean-curvature flow step.
 */

#ifndef GEOMLIB_FAIRING_H
#define GEOMLIB_FAIRING_H

#include "geomlib_mesh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float lambda;          /* Flow time step (default 0.1) */
    int iterations;        /* Number of flow steps (default 5) */
    bool preserve_boundaries;
} geomlib_fairing_params;

geomlib_fairing_params geomlib_fairing_params_default(void);

/* Smooth vertex positions in place with implicit mean-curvature flow. */
int geomlib_fair_mesh(geomlib_trimesh *mesh,
                      const geomlib_fairing_params *params);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_FAIRING_H */
