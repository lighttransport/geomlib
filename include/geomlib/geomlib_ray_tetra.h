/*
 * geomlib - Ray / tetrahedron intersection
 *
 * Returns the nearest triangle hit among the four tetrahedron faces.
 */

#ifndef GEOMLIB_RAY_TETRA_H
#define GEOMLIB_RAY_TETRA_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;
    float t;             /* distance along `direction` */
    geomlib_vec3 point;
    geomlib_vec3 normal; /* geometric face normal (follows face winding) */
} geomlib_ray_tetra_hit;

int geomlib_ray_tetrahedron_intersect(
    geomlib_vec3 origin, geomlib_vec3 direction,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2, geomlib_vec3 v3,
    geomlib_ray_tetra_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RAY_TETRA_H */
