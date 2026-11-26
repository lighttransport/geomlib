/*
 * geomlib - Linear blend skinning (see geomlib_skinning.h)
 *
 * Clean-room implementation of standard linear blend skinning (LBS):
 *   x' = sum_i w_i * T_i * x, normalized by sum_i w_i.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_skinning.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

int geomlib_skin_mesh(const geomlib_trimesh *input,
                      const geomlib_skin_weight *weights,
                      const geomlib_mat4 *bone_transforms,
                      size_t bone_count,
                      geomlib_trimesh *output)
{
    geomlib_trimesh result = {0};
    size_t i, k;
    int rc;

    if (!input || !output || !input->vertices || !input->indices ||
        !weights || !bone_transforms || bone_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    rc = geomlib_trimesh_init(&result, input->vertex_count, input->triangle_count);
    if (rc != GEOMLIB_MESH_SUCCESS) return rc;

    for (i = 0; i < input->vertex_count; ++i) {
        double x = 0.0, y = 0.0, z = 0.0, wsum = 0.0;
        for (k = 0; k < GEOMLIB_SKIN_MAX_BONES; ++k) {
            float w = weights[i].weights[k];
            uint32_t b = weights[i].bone_indices[k];
            geomlib_vec4 p;
            if (!(w > 0.0f) || b >= bone_count) continue;
            p = geomlib_mat4_mul_vec4(bone_transforms[b],
                geomlib_vec4_create(input->vertices[i].x, input->vertices[i].y,
                                    input->vertices[i].z, 1.0f));
            x += w*p.x;
            y += w*p.y;
            z += w*p.z;
            wsum += w;
        }
        if (wsum > 1e-12) {
            geomlib_trimesh_add_vertex(&result,
                geomlib_vec3_create((float)(x/wsum), (float)(y/wsum), (float)(z/wsum)),
                NULL, NULL);
        } else {
            geomlib_trimesh_add_vertex(&result, input->vertices[i], NULL, NULL);
        }
    }
    for (i = 0; i < input->triangle_count; ++i)
        geomlib_trimesh_add_triangle(&result, input->indices[i*3],
                                     input->indices[i*3+1], input->indices[i*3+2]);
    *output = result;
    return GEOMLIB_MESH_SUCCESS;
}
