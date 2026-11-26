/*
 * geomlib - Ray / infinite cone intersection
 *
 * Analytic intersection of a ray with one nappe of a right circular cone.
 */

#ifndef GEOMLIB_RAY_CONE_H
#define GEOMLIB_RAY_CONE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;
    double t;             /* parameter along `direction` */
    geomlib_vec3 point;   /* surface hit point */
} geomlib_ray_cone_hit;

/* Cone: apex, unit axis direction, half-angle in radians (0 < angle < pi/2).
 * Only the nappe in the axis direction is considered. `direction` need not be
 * normalized. Returns the nearest hit with t >= 0. */
int geomlib_ray_cone_intersect(geomlib_vec3 origin, geomlib_vec3 direction,
                               geomlib_vec3 apex, geomlib_vec3 axis,
                               float half_angle, geomlib_ray_cone_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RAY_CONE_H */
