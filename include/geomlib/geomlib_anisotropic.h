/*
 * geomlib - Anisotropic curvature-adaptive remeshing
 *
 * Uses a per-vertex symmetric metric tensor (2x2 in the tangent plane) to
 * measure edge length; remeshing splits/collapses/flips edges against metric
 * length targets so vertices concentrate in high-curvature regions.
 */

#ifndef GEOMLIB_ANISOTROPIC_H
#define GEOMLIB_ANISOTROPIC_H

#include "geomlib_mesh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float target_length;       /* Base metric edge length */
    float curvature_scale;     /* How strongly curvature refines (default 2.0) */
    int iterations;            /* Remeshing iterations (default 5) */
    bool preserve_features;    /* Freeze sharp edges (default true) */
    float feature_angle;       /* Dihedral angle threshold in degrees */
} geomlib_anisotropic_params;

geomlib_anisotropic_params geomlib_anisotropic_params_default(float target_length);

/*
 * Remesh in place using curvature-derived anisotropic metric tensors.
 * The output keeps the input topology style (triangle mesh).
 */
int geomlib_remesh_anisotropic(geomlib_trimesh *mesh,
                               const geomlib_anisotropic_params *params);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ANISOTROPIC_H */
