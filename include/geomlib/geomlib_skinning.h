/*
 * geomlib - Linear blend skinning
 *
 * Deforms mesh vertices by a weighted blend of rigid bone transforms.
 */

#ifndef GEOMLIB_SKINNING_H
#define GEOMLIB_SKINNING_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GEOMLIB_SKIN_MAX_BONES 4

typedef struct {
    uint32_t bone_indices[GEOMLIB_SKIN_MAX_BONES];
    float weights[GEOMLIB_SKIN_MAX_BONES];
} geomlib_skin_weight;

/* Apply bone transforms to a mesh using per-vertex skin weights. Each bone
 * transform maps from rest space to deformed space. */
int geomlib_skin_mesh(const geomlib_trimesh *input,
                      const geomlib_skin_weight *weights,
                      const geomlib_mat4 *bone_transforms,
                      size_t bone_count,
                      geomlib_trimesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SKINNING_H */
