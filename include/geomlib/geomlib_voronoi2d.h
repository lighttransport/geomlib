/*
 * geomlib - 2D Voronoi diagram (bounded by a user-supplied box)
 *
 * Each cell is the intersection of the bounding box with the perpendicular
 * bisector half-planes between its site and every other site.
 */

#ifndef GEOMLIB_VORONOI2D_H
#define GEOMLIB_VORONOI2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *sites;      /* owned copy of de-duplicated input points */
    size_t site_count;
    geomlib_vec2 *cells;      /* concatenated convex polygons (CCW) */
    size_t *offsets;          /* site_count+1 offsets into `cells` */
    size_t vertex_count;
    geomlib_vec2 bounds_min, bounds_max;
} geomlib_voronoi2d;

geomlib_voronoi2d geomlib_voronoi2d_create(void);
void geomlib_voronoi2d_free(geomlib_voronoi2d *v);

/* Build clipped Voronoi cells. Duplicate sites (within 1e-9) are ignored. */
int geomlib_compute_voronoi2d(const geomlib_vec2 *points, size_t count,
                              geomlib_vec2 bounds_min, geomlib_vec2 bounds_max,
                              geomlib_voronoi2d *out);

/* Index of the site whose cell contains q, or SIZE_MAX when outside the box.
 * For a clipped diagram this also equals the nearest-site query. */
size_t geomlib_voronoi2d_locate(const geomlib_voronoi2d *v, geomlib_vec2 q);

/* Area of one cell (0 for invalid index). */
double geomlib_voronoi2d_cell_area(const geomlib_voronoi2d *v, size_t i);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_VORONOI2D_H */
