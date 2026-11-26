/*
 * geomlib - 2D polygon offsetting
 *
 * Offsets a simple polygon inward or outward by moving each edge along its
 * normal and intersecting consecutive offset lines. Reflex vertices are
 * handled with round joins (arc approximation); self-intersections of the
 * offset are removed by a deterministic loop-trimming pass.
 */

#ifndef GEOMLIB_POLYGON_OFFSET_H
#define GEOMLIB_POLYGON_OFFSET_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_polygon2d *polygons;  /* Result loops (may be multiple) */
    size_t count;
} geomlib_offset_result;

geomlib_offset_result geomlib_offset_result_create(void);
void geomlib_offset_result_free(geomlib_offset_result *r);

/*
 * Offset a simple 2D polygon. Positive distance offsets outward, negative
 * inward. The input is given in CCW order.
 */
int geomlib_offset_polygon(const geomlib_polygon2d *polygon,
                           float distance,
                           int arc_segments,
                           geomlib_offset_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_OFFSET_H */
