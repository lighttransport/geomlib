/*
 * geomlib - Alpha wrapping
 *
 * Produces a watertight, intersection-free enclosing surface from arbitrary
 * input geometry by offsetting a signed-distance representation and extracting
 * the offset isosurface.
 */

#ifndef GEOMLIB_ALPHA_WRAP_H
#define GEOMLIB_ALPHA_WRAP_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t resolution;   /* SDF grid resolution (default 48) */
    float alpha;           /* Offset distance (default 0.01 * bbox diagonal) */
    int smoothing_steps;   /* Optional Gaussian smoothing (default 1) */
} geomlib_alpha_wrap_params;

geomlib_alpha_wrap_params geomlib_alpha_wrap_params_default(void);

/*
 * Wrap `input` with a watertight surface. `input` may be open or contain
 * duplicated faces; the SDF representation handles it robustly.
 */
int geomlib_alpha_wrap(const geomlib_trimesh *input,
                       const geomlib_alpha_wrap_params *params,
                       geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ALPHA_WRAP_H */
