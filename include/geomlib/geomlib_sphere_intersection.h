/*
 * geomlib - Sphere-sphere intersection
 *
 * Computes the circle where two spheres intersect (or detects containment,
 * tangency and disjointness).
 */

#ifndef GEOMLIB_SPHERE_INTERSECTION_H
#define GEOMLIB_SPHERE_INTERSECTION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_SPHERE_NONE = 0,
    GEOMLIB_SPHERE_POINT,
    GEOMLIB_SPHERE_CIRCLE,
    GEOMLIB_SPHERE_CONTAINED
} geomlib_sphere_intersect_kind;

typedef struct {
    geomlib_sphere_intersect_kind kind;
    geomlib_vec3 center;   /* centre of the intersection circle / point */
    geomlib_vec3 normal;   /* circle axis (unit) */
    float radius;
} geomlib_sphere_intersection;

geomlib_sphere_intersection geomlib_sphere_sphere_intersection(
    geomlib_vec3 c0, float r0, geomlib_vec3 c1, float r1);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPHERE_INTERSECTION_H */
