/*
 * geomlib - Distance between two triangles
 *
 * Returns the minimum distance and closest points, using vertex-to-triangle
 * projections and edge-to-edge segment distances.
 */

#ifndef GEOMLIB_TRI_TRI_DISTANCE_H
#define GEOMLIB_TRI_TRI_DISTANCE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float distance;         /* 0 when the triangles intersect */
    geomlib_vec3 point_a;   /* closest point on triangle A */
    geomlib_vec3 point_b;   /* closest point on triangle B */
} geomlib_tri_tri_distance_result;

int geomlib_triangle_triangle_distance(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2,
    geomlib_vec3 b0, geomlib_vec3 b1, geomlib_vec3 b2,
    geomlib_tri_tri_distance_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TRI_TRI_DISTANCE_H */
