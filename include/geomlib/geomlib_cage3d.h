/*
 * geomlib - Tetrahedral cage deformation
 *
 * Deforms points with affine-precise barycentric coordinates of a tetrahedron.
 */

#ifndef GEOMLIB_CAGE3D_H
#define GEOMLIB_CAGE3D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* `cage` and `deformed` are 4 vertices each. `points`/`out` have `count`
 * entries. Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_deform_tetra_cage(const geomlib_vec3 cage[4],
                              const geomlib_vec3 deformed[4],
                              const geomlib_vec3 *points, size_t count,
                              geomlib_vec3 *out);

/* Barycentric coordinates of q in the cage (4 weights summing to 1). */
int geomlib_tetra_weights(const geomlib_vec3 cage[4], geomlib_vec3 q,
                          double weights[4]);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CAGE3D_H */
