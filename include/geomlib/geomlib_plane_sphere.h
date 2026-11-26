/*
 * geomlib - Plane / sphere intersection
 *
 * Computes the circle where a plane cuts a sphere, with tangency and disjoint
 * classification.
 */

#ifndef GEOMLIB_PLANE_SPHERE_H
#define GEOMLIB_PLANE_SPHERE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_PLANE_SPHERE_NONE = 0,
    GEOMLIB_PLANE_SPHERE_TANGENT,
    GEOMLIB_PLANE_SPHERE_CIRCLE
} geomlib_plane_sphere_kind;

typedef struct {
    geomlib_plane_sphere_kind kind;
    geomlib_vec3 center;   /* centre of the intersection circle */
    geomlib_vec3 normal;   /* unit plane normal */
    float radius;          /* circle radius */
} geomlib_plane_sphere_intersection;

geomlib_plane_sphere_intersection geomlib_plane_sphere_intersect(
    geomlib_vec3 plane_point, geomlib_vec3 plane_normal,
    geomlib_vec3 sphere_center, float radius);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PLANE_SPHERE_H */
