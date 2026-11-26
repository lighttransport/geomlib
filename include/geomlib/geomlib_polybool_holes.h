/*
 * geomlib - 2D polygon Boolean operations with hole support
 *
 * Bridges holes into simple polygons, delegates to the existing Boolean
 * engine, then removes degenerate bridge fragments.
 */

#ifndef GEOMLIB_POLYBOOL_HOLES_H
#define GEOMLIB_POLYBOOL_HOLES_H

#include "geomlib_polybool2d.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *vertices;   /* owned: outer then holes */
    size_t *offsets;          /* hole_count+1 offsets into vertices */
    size_t hole_count;
} geomlib_polygon_with_holes;

geomlib_polygon_with_holes geomlib_polygon_with_holes_create(void);
void geomlib_polygon_with_holes_free(geomlib_polygon_with_holes *p);

/* Build from an outer ring and `nh` hole rings (offsets computed internally). */
int geomlib_polygon_with_holes_set(const geomlib_vec2 *outer, size_t no,
                                   const geomlib_vec2 *const *holes,
                                   const size_t *hole_sizes, size_t nh,
                                   geomlib_polygon_with_holes *out);

/* Boolean operation on two polygons with holes. */
int geomlib_polygon_boolean_holes(const geomlib_polygon_with_holes *a,
                                  const geomlib_polygon_with_holes *b,
                                  geomlib_poly_bool_op op,
                                  geomlib_poly_bool_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYBOOL_HOLES_H */
