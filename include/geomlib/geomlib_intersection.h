/*
 * geomlib - Triangle-mesh intersection
 *
 * Computes the intersection polyline of two triangle meshes as a set of
 * segments. Each segment joins two crossing points on the same pair of
 * faces (one from each mesh).
 */

#ifndef GEOMLIB_INTERSECTION_H
#define GEOMLIB_INTERSECTION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *segments;   /* two endpoints per segment */
    size_t segment_count;
} geomlib_intersection;

geomlib_intersection geomlib_intersection_create(void);
void geomlib_intersection_free(geomlib_intersection *r);

/* Intersect mesh `a` against mesh `b`. Coplanar overlapping faces are not
 * reported. Returns an empty result on allocation failure. */
geomlib_intersection geomlib_mesh_mesh_intersection(const geomlib_trimesh *a,
                                                     const geomlib_trimesh *b);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_INTERSECTION_H */
