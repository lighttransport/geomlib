/*
 * geomlib - Sphere / capsule overlap
 *
 * Compares the sphere-centre-to-axis distance with the sum of radii.
 */

#ifndef GEOMLIB_SPHERE_CAPSULE_H
#define GEOMLIB_SPHERE_CAPSULE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;       /* 1 on touch or overlap */
    float distance;       /* axis distance minus the sum of radii */
    geomlib_vec3 point;   /* closest point on the capsule axis */
} geomlib_sphere_capsule_result;

int geomlib_sphere_capsule_overlap(geomlib_vec3 center, float radius,
                                   geomlib_vec3 a, geomlib_vec3 b,
                                   float cap_radius,
                                   geomlib_sphere_capsule_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPHERE_CAPSULE_H */
