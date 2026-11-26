/*
 * geomlib - Capsule / oriented bounding box overlap
 *
 * Transforms the capsule axis into the OBB frame and reuses the segment-box
 * distance.
 */

#ifndef GEOMLIB_CAPSULE_OBB_H
#define GEOMLIB_CAPSULE_OBB_H

#include "geomlib_mesh.h"
#include "geomlib_collision.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;   /* 1 on touch or overlap */
    float distance;   /* axis-to-OBB distance minus the capsule radius */
} geomlib_capsule_obb_result;

int geomlib_capsule_obb_overlap(geomlib_vec3 a, geomlib_vec3 b, float r,
                                const geomlib_obb *obb,
                                geomlib_capsule_obb_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CAPSULE_OBB_H */
