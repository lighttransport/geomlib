/*
 * geomlib - Capsule / triangle overlap
 *
 * Uses the exact distance between the capsule axis segment and the triangle.
 */

#ifndef GEOMLIB_CAPSULE_TRIANGLE_H
#define GEOMLIB_CAPSULE_TRIANGLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float distance;          /* minimum distance */
    geomlib_vec3 point_a;    /* closest point on the segment */
    geomlib_vec3 point_b;    /* closest point on the triangle */
} geomlib_segment_triangle_distance_result;

int geomlib_segment_triangle_distance(
    geomlib_vec3 a, geomlib_vec3 b,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_segment_triangle_distance_result *out);

typedef struct {
    int intersects;          /* 1 on touch or overlap */
    float distance;          /* axis-to-triangle distance minus radius */
    geomlib_vec3 point_a;    /* closest point on the capsule axis */
    geomlib_vec3 point_b;    /* closest point on the triangle */
} geomlib_capsule_triangle_result;

int geomlib_capsule_triangle_overlap(
    geomlib_vec3 a, geomlib_vec3 b, float r,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_capsule_triangle_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CAPSULE_TRIANGLE_H */
