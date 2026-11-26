/*
 * geomlib - Rigid/similarity point-set fitting (Kabsch/Horn)
 *
 * Fits a rotation, translation and optional uniform scale between two
 * corresponding point sets.
 */

#ifndef GEOMLIB_RIGID_FIT_H
#define GEOMLIB_RIGID_FIT_H

#include "geomlib_mesh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Fit src -> dst. `transform` is column-major. Returns GEOMLIB_MESH_SUCCESS. */
int geomlib_fit_rigid_points(const geomlib_vec3 *src, const geomlib_vec3 *dst,
                             size_t count, bool allow_scale,
                             geomlib_mat4 *transform, float *rms_out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RIGID_FIT_H */
