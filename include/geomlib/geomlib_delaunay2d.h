/*
 * geomlib - 2D Delaunay triangulation
 *
 * Incremental Bowyer-Watson construction for 2D point sets, returning the
 * triangle index list. Uses exact-ish double predicates.
 */

#ifndef GEOMLIB_DELAUNAY2D_H
#define GEOMLIB_DELAUNAY2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *triangles;   /* 3 indices per triangle */
    size_t triangle_count;
} geomlib_delaunay2d_result;

geomlib_delaunay2d_result geomlib_delaunay2d_result_create(void);
void geomlib_delaunay2d_result_free(geomlib_delaunay2d_result *r);

/* Triangulate a set of 2D points (duplicates are ignored). */
int geomlib_delaunay2d(const geomlib_vec2 *points, size_t count,
                       geomlib_delaunay2d_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DELAUNAY2D_H */
