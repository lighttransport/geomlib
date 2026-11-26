/*
 * geomlib - Capsule collision primitives
 *
 * A capsule is a line segment with a radius; its surface is the Minkowski sum
 * of the segment with a sphere.
 */

#ifndef GEOMLIB_CAPSULE_H
#define GEOMLIB_CAPSULE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;          /* 1 when the capsules touch or overlap */
    float distance;          /* axis distance minus the sum of radii */
    geomlib_vec3 point_a;    /* closest point on axis A */
    geomlib_vec3 point_b;    /* closest point on axis B */
} geomlib_capsule_result;

/* Overlap test for capsules (a0,a1,r0) and (b0,b1,r1). */
int geomlib_capsule_capsule_overlap(geomlib_vec3 a0, geomlib_vec3 a1, float r0,
                                    geomlib_vec3 b0, geomlib_vec3 b1, float r1,
                                    geomlib_capsule_result *out);

/* Signed distance from `q` to the surface of capsule (a0,a1,r): negative when
 * inside. The closest point on the axis is written to `axis_point` if given. */
float geomlib_capsule_point_distance(geomlib_vec3 q, geomlib_vec3 a0,
                                     geomlib_vec3 a1, float r,
                                     geomlib_vec3 *axis_point);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CAPSULE_H */
