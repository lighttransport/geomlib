/*
 * geomlib - Tetrahedron / tetrahedron intersection
 *
 * Closest-feature test for two tetrahedra.
 */

#ifndef GEOMLIB_TETRA_TETRA_H
#define GEOMLIB_TETRA_TETRA_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

int geomlib_tetrahedra_intersect(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2, geomlib_vec3 a3,
    geomlib_vec3 b0, geomlib_vec3 b1, geomlib_vec3 b2, geomlib_vec3 b3);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TETRA_TETRA_H */
