/*
 * geomlib - 2D convex hull
 *
 * Andrew's monotone chain: sort by (x,y), build lower and upper chains with
 * cross-product popping, producing a CCW hull without duplicate points.
 */

#ifndef GEOMLIB_HULL2D_H
#define GEOMLIB_HULL2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *vertices;   /* CCW hull vertices, no repeated closure point */
    size_t count;
} geomlib_hull2d;

geomlib_hull2d geomlib_hull2d_create(void);
void geomlib_hull2d_free(geomlib_hull2d *h);

/* Convex hull of 2D points. Duplicates and collinear points are removed. */
int geomlib_convex_hull2d(const geomlib_vec2 *points, size_t count,
                          geomlib_hull2d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_HULL2D_H */
