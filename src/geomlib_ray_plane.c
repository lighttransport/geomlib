/*
 * geomlib - Ray/plane intersection (see geomlib_ray_plane.h)
 *
 * Clean-room plane equation solve. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_ray_plane.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_ray_plane_intersection(geomlib_vec3 origin, geomlib_vec3 direction,
                                   geomlib_vec3 plane_point,
                                   geomlib_vec3 plane_normal,
                                   float *t_out, geomlib_vec3 *point_out)
{
    float den=geomlib_vec3_dot(plane_normal,direction);
    float t;
    if (fabsf(den)<1e-9f) return 0;
    t=geomlib_vec3_dot(plane_normal,geomlib_vec3_sub(plane_point,origin))/den;
    if (t<0.0f) return 0;
    if (t_out) *t_out=t;
    if (point_out) *point_out=geomlib_vec3_add(origin,geomlib_vec3_mul(direction,t));
    return 1;
}

int geomlib_segment_plane_intersection(geomlib_vec3 a, geomlib_vec3 b,
                                       geomlib_vec3 plane_point,
                                       geomlib_vec3 plane_normal,
                                       float *t_out, geomlib_vec3 *point_out)
{
    float den=geomlib_vec3_dot(plane_normal,geomlib_vec3_sub(b,a));
    float t;
    if (fabsf(den)<1e-9f) return 0;
    t=geomlib_vec3_dot(plane_normal,geomlib_vec3_sub(plane_point,a))/den;
    if (t<0.0f || t>1.0f) return 0;
    if (t_out) *t_out=t;
    if (point_out) *point_out=geomlib_vec3_add(a,geomlib_vec3_mul(geomlib_vec3_sub(b,a),t));
    return 1;
}
