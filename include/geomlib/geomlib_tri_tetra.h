/*
 * geomlib - Triangle / tetrahedron intersection
 *
 * Reports whether a triangle and a tetrahedron have any common point.
 */

#ifndef GEOMLIB_TRI_TETRA_H
#define GEOMLIB_TRI_TETRA_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

int geomlib_triangle_tetrahedron_intersect(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2, geomlib_vec3 v3);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TRI_TETRA_H */
