/*
 * geomlib - Ray/plane and segment/plane intersection
 *
 * Computes the intersection point of a ray with an infinite plane.
 */

#ifndef GEOMLIB_RAY_PLANE_H
#define GEOMLIB_RAY_PLANE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Returns 1 when the ray meets the plane at t >= 0, else 0. */
int geomlib_ray_plane_intersection(geomlib_vec3 origin, geomlib_vec3 direction,
                                   geomlib_vec3 plane_point,
                                   geomlib_vec3 plane_normal,
                                   float *t_out, geomlib_vec3 *point_out);

/* Returns 1 when the segment meets the plane (t in [0,1]). */
int geomlib_segment_plane_intersection(geomlib_vec3 a, geomlib_vec3 b,
                                       geomlib_vec3 plane_point,
                                       geomlib_vec3 plane_normal,
                                       float *t_out, geomlib_vec3 *point_out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RAY_PLANE_H */
