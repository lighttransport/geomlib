/*
 * geomlib - Segment / triangle intersection
 *
 * Möller-Trumbore intersection of a finite 3D segment with a triangle.
 */

#ifndef GEOMLIB_SEGMENT_TRIANGLE_H
#define GEOMLIB_SEGMENT_TRIANGLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;
    float t;             /* parameter along the segment, in [0,1] */
    geomlib_vec3 point;  /* intersection point */
    geomlib_vec3 bary;   /* (b0,b1,b2) barycentric coordinates */
} geomlib_segment_triangle_hit;

/* Returns 1 when the segment meets the triangle. Coplanar segments are not
 * reported (the intersection is not a single point). */
int geomlib_segment_triangle_intersect(
    geomlib_vec3 a, geomlib_vec3 b,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2,
    geomlib_segment_triangle_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SEGMENT_TRIANGLE_H */
