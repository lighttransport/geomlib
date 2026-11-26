/*
 * geomlib - Sphere / AABB overlap
 *
 * Clamps the sphere centre to the box and compares the squared distance with
 * the squared radius.
 */

#ifndef GEOMLIB_SPHERE_AABB_H
#define GEOMLIB_SPHERE_AABB_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;          /* 1 on overlap or tangency */
    geomlib_vec3 closest;    /* closest point on the box to the sphere centre */
    float distance;          /* centre-to-box distance minus radius (<0 when
                              * the sphere penetrates the box) */
} geomlib_sphere_aabb_result;

/* Returns 1 when a sphere overlaps the axis-aligned box [bmin,bmax]. */
int geomlib_sphere_aabb_overlap(geomlib_vec3 center, float radius,
                                geomlib_vec3 bmin, geomlib_vec3 bmax,
                                geomlib_sphere_aabb_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPHERE_AABB_H */
