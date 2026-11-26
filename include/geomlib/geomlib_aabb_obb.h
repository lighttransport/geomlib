/*
 * geomlib - AABB / oriented bounding box overlap
 *
 * Separating-axis test over 3 box axes, 3 OBB axes and 9 cross-product axes.
 */

#ifndef GEOMLIB_AABB_OBB_H
#define GEOMLIB_AABB_OBB_H

#include "geomlib_mesh.h"
#include "geomlib_collision.h"

#ifdef __cplusplus
extern "C" {
#endif

int geomlib_aabb_obb_overlap(geomlib_vec3 amin, geomlib_vec3 amax,
                             const geomlib_obb *obb);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_AABB_OBB_H */
