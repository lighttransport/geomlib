/*
 * geomlib - Polyline / triangle-mesh intersection
 *
 * Traces every segment of a polyline through a BVH and returns all crossing
 * points on the mesh surface.
 */

#ifndef GEOMLIB_POLYLINE_INTERSECTION_H
#define GEOMLIB_POLYLINE_INTERSECTION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *points;   /* crossing points along the polyline */
    size_t count;
} geomlib_polyline_hits;

geomlib_polyline_hits geomlib_polyline_hits_create(void);
void geomlib_polyline_hits_free(geomlib_polyline_hits *h);

int geomlib_polyline_mesh_intersection(const geomlib_vec3 *polyline, size_t n,
                                       const geomlib_trimesh *mesh,
                                       geomlib_polyline_hits *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYLINE_INTERSECTION_H */
