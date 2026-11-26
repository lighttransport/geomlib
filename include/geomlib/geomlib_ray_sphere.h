/*
 * geomlib - Ray / sphere intersection
 *
 * Analytic quadratic roots along a ray for a sphere in 3D.
 */

#ifndef GEOMLIB_RAY_SPHERE_H
#define GEOMLIB_RAY_SPHERE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double t0, t1;             /* near and far roots (t0 <= t1) */
    geomlib_vec3 p0, p1;       /* entry / exit points */
    geomlib_vec3 n0, n1;       /* outward surface normals there */
} geomlib_ray_sphere_hit;

/* Returns 1 when the ray meets the sphere, else 0. `direction` need not be
 * normalized; roots are measured in units of `direction`. */
int geomlib_ray_sphere_intersect(geomlib_vec3 origin, geomlib_vec3 direction,
                                 geomlib_vec3 center, float radius,
                                 geomlib_ray_sphere_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RAY_SPHERE_H */
