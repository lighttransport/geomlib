/*
 * geomlib - Distance between two axis-aligned boxes
 *
 * Returns the minimum Euclidean distance and closest points on each box.
 */

#ifndef GEOMLIB_AABB_DISTANCE_H
#define GEOMLIB_AABB_DISTANCE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float distance;          /* 0 when the boxes overlap or touch */
    geomlib_vec3 point_a;    /* closest point on box A */
    geomlib_vec3 point_b;    /* closest point on box B */
} geomlib_aabb_distance_result;

int geomlib_aabb_distance(geomlib_vec3 amin, geomlib_vec3 amax,
                          geomlib_vec3 bmin, geomlib_vec3 bmax,
                          geomlib_aabb_distance_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_AABB_DISTANCE_H */
