/*
 * geomlib - As-rigid-as-possible surface deformation
 *
 * Local/global ARAP deformation of a triangle mesh with positional handles.
 */

#ifndef GEOMLIB_DEFORM_H
#define GEOMLIB_DEFORM_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t vertex;        /* Handle vertex index */
    geomlib_vec3 target;    /* Desired deformed position */
} geomlib_deform_handle;

/*
 * Deform `input` so its handle vertices move to their targets while the rest
 * of the surface deforms as rigidly as possible. `output` is a new mesh with
 * the same topology and the deformed vertex positions.
 */
int geomlib_deform_arap(const geomlib_trimesh *input,
                        const geomlib_deform_handle *handles,
                        size_t handle_count,
                        int iterations,
                        geomlib_trimesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DEFORM_H */
