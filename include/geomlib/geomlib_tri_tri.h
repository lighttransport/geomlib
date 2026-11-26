/*
 * geomlib - Triangle / triangle intersection
 *
 * Returns the intersection segment of two triangles, including coplanar cases.
 */

#ifndef GEOMLIB_TRI_TRI_H
#define GEOMLIB_TRI_TRI_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;        /* 1 when the triangles intersect */
    geomlib_vec3 p0, p1;   /* intersection segment endpoints */
} geomlib_tri_tri_result;

/* Intersection of triangles (a0,a1,a2) and (b0,b1,b2). Winding irrelevant. */
geomlib_tri_tri_result geomlib_triangle_intersection(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2,
    geomlib_vec3 b0, geomlib_vec3 b1, geomlib_vec3 b2);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TRI_TRI_H */
