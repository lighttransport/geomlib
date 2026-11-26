/*
 * geomlib - Closest points between two 3D line segments
 *
 * Computes the minimum distance between two finite segments, returning the
 * closest point on each segment and the segment parameters.
 */

#ifndef GEOMLIB_SEGMENT_DISTANCE_H
#define GEOMLIB_SEGMENT_DISTANCE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 point_a, point_b;
    float s, t;        /* parameters in [0,1] along each segment */
    float distance;
} geomlib_segment_distance;

int geomlib_segment_closest_points(geomlib_vec3 a0, geomlib_vec3 a1,
                                   geomlib_vec3 b0, geomlib_vec3 b1,
                                   geomlib_segment_distance *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SEGMENT_DISTANCE_H */
