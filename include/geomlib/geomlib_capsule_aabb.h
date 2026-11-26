/*
 * geomlib - Capsule / AABB overlap
 *
 * Uses the exact minimum distance between the capsule axis segment and the box.
 */

#ifndef GEOMLIB_CAPSULE_AABB_H
#define GEOMLIB_CAPSULE_AABB_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Minimum Euclidean distance from segment [a,b] to the box [bmin,bmax]. */
float geomlib_segment_aabb_distance(geomlib_vec3 a, geomlib_vec3 b,
                                    geomlib_vec3 bmin, geomlib_vec3 bmax);

typedef struct {
    int intersects;   /* 1 when the capsule touches or overlaps the box */
    float distance;   /* axis-to-box distance minus radius */
} geomlib_capsule_aabb_result;

int geomlib_capsule_aabb_overlap(geomlib_vec3 a, geomlib_vec3 b, float r,
                                 geomlib_vec3 bmin, geomlib_vec3 bmax,
                                 geomlib_capsule_aabb_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CAPSULE_AABB_H */
