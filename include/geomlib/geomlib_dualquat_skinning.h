/*
 * geomlib - Dual-quaternion skinning
 *
 * Blends rigid transforms as dual quaternions and normalizes the result,
 * avoiding the shrinking artefacts of linear blend skinning.
 */

#ifndef GEOMLIB_DUALQUAT_SKINNING_H
#define GEOMLIB_DUALQUAT_SKINNING_H

#include "geomlib_mesh.h"
#include "geomlib_skinning.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Skin a mesh with dual-quaternion blending of bone transforms. */
int geomlib_skin_mesh_dualquat(const geomlib_trimesh *input,
                               const geomlib_skin_weight *weights,
                               const geomlib_mat4 *bone_transforms,
                               size_t bone_count,
                               geomlib_trimesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DUALQUAT_SKINNING_H */
