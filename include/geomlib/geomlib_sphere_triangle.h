/*
 * geomlib - Sphere / triangle overlap
 *
 * Tests whether a sphere touches or penetrates a triangle by comparing the
 * radius with the distance from the centre to the closest triangle point.
 */

#ifndef GEOMLIB_SPHERE_TRIANGLE_H
#define GEOMLIB_SPHERE_TRIANGLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;         /* 1 on overlap or tangency */
    geomlib_vec3 closest;   /* closest point on the triangle */
    float distance;         /* centre-to-triangle distance minus radius */
} geomlib_sphere_triangle_result;

int geomlib_sphere_triangle_overlap(
    geomlib_vec3 center, float radius,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_sphere_triangle_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPHERE_TRIANGLE_H */
