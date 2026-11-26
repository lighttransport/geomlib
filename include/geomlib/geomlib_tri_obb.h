/*
 * geomlib - Triangle / oriented bounding box overlap
 *
 * Transforms the triangle into the OBB frame and reuses the triangle-AABB SAT.
 */

#ifndef GEOMLIB_TRI_OBB_H
#define GEOMLIB_TRI_OBB_H

#include "geomlib_mesh.h"
#include "geomlib_collision.h"

#ifdef __cplusplus
extern "C" {
#endif

int geomlib_triangle_obb_overlap(geomlib_vec3 a0, geomlib_vec3 a1,
                                 geomlib_vec3 a2, const geomlib_obb *obb);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TRI_OBB_H */
