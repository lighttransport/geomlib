/*
 * geomlib - Ray / oriented bounding box intersection
 *
 * Slab test performed in the OBB's local orthonormal frame.
 */

#ifndef GEOMLIB_RAY_OBB_H
#define GEOMLIB_RAY_OBB_H

#include "geomlib_mesh.h"
#include "geomlib_collision.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double t;             /* distance parameter along `direction` */
    geomlib_vec3 point;   /* surface hit point */
    geomlib_vec3 normal;  /* outward face normal in world space */
} geomlib_ray_obb_hit;

/* Returns 1 for the nearest hit with t >= 0. */
int geomlib_ray_obb_intersect(geomlib_vec3 origin, geomlib_vec3 direction,
                              const geomlib_obb *obb,
                              geomlib_ray_obb_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RAY_OBB_H */
