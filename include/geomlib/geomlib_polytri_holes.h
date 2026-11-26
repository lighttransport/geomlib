/*
 * geomlib - 2D polygon triangulation with holes
 *
 * Bridges holes into the outer boundary, then ear-clips the resulting simple
 * polygon. Triangle indices refer to concatenated [outer..., hole0..., ...].
 */

#ifndef GEOMLIB_POLYTRI_HOLES_H
#define GEOMLIB_POLYTRI_HOLES_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *triangles;   /* 3 indices per triangle */
    size_t triangle_count;
    geomlib_vec2 *vertices; /* owned concatenation of outer and hole points */
    size_t vertex_count;
} geomlib_polytri_result;

geomlib_polytri_result geomlib_polytri_result_create(void);
void geomlib_polytri_result_free(geomlib_polytri_result *r);

/* Triangulate an outer polygon with holes. `holes` is an array of `nh` arrays
 * of `hole_sizes[i]` points; winding is normalized internally. */
int geomlib_triangulate_with_holes(const geomlib_vec2 *outer, size_t no,
                                   const geomlib_vec2 *const *holes,
                                   const size_t *hole_sizes, size_t nh,
                                   geomlib_polytri_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYTRI_HOLES_H */
