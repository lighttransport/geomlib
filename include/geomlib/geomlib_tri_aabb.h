/*
 * geomlib - Triangle / AABB overlap
 *
 * Separating-axis test between a triangle and an axis-aligned box.
 */

#ifndef GEOMLIB_TRI_AABB_H
#define GEOMLIB_TRI_AABB_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Returns 1 when triangle (v0,v1,v2) overlaps the box [bmin,bmax]. */
int geomlib_triangle_aabb_overlap(geomlib_vec3 v0, geomlib_vec3 v1,
                                  geomlib_vec3 v2,
                                  geomlib_vec3 bmin, geomlib_vec3 bmax);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TRI_AABB_H */
