/*
 * geomlib - Sphere / oriented bounding box overlap
 *
 * Clamps the sphere centre to the OBB in its local frame.
 */

#ifndef GEOMLIB_SPHERE_OBB_H
#define GEOMLIB_SPHERE_OBB_H

#include "geomlib_mesh.h"
#include "geomlib_collision.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;         /* 1 on overlap or tangency */
    float distance;         /* centre-to-box distance minus radius */
    geomlib_vec3 closest;   /* closest point on the OBB to the sphere centre */
} geomlib_sphere_obb_result;

int geomlib_sphere_obb_overlap(geomlib_vec3 center, float radius,
                               const geomlib_obb *obb,
                               geomlib_sphere_obb_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPHERE_OBB_H */
