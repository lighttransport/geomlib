/*
 * geomlib - Ray / infinite cylinder intersection
 *
 * Analytic intersection of a ray with a cylinder defined by an axis and radius.
 */

#ifndef GEOMLIB_RAY_CYLINDER_H
#define GEOMLIB_RAY_CYLINDER_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double t;             /* distance parameter along `direction` */
    geomlib_vec3 point;   /* surface hit point */
    geomlib_vec3 normal;  /* outward surface normal at the hit */
} geomlib_ray_cylinder_hit;

/* Returns 1 for the nearest hit with t >= 0. Rays parallel to the axis are not
 * reported. `direction` need not be normalized. */
int geomlib_ray_cylinder_intersect(
    geomlib_vec3 origin, geomlib_vec3 direction,
    geomlib_vec3 axis_point, geomlib_vec3 axis_dir, float radius,
    geomlib_ray_cylinder_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RAY_CYLINDER_H */
