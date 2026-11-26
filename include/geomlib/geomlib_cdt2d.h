/*
 * geomlib - 2D constrained Delaunay triangulation
 *
 * Builds a Delaunay triangulation that is forced to contain a set of input
 * segments. Segments are recovered with edge flips and, when a quadrilateral
 * is not flippable, by inserting a Steiner point at the crossing.
 */

#ifndef GEOMLIB_CDT2D_H
#define GEOMLIB_CDT2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *vertices;      /* owned: input points plus Steiner points */
    size_t vertex_count;
    uint32_t *triangles;         /* 3 indices per triangle */
    size_t triangle_count;
    uint32_t *constraints;       /* 2 indices per recovered constraint */
    size_t constraint_count;
} geomlib_cdt2d;

geomlib_cdt2d geomlib_cdt2d_create(void);
void geomlib_cdt2d_free(geomlib_cdt2d *c);

/* Triangulate `points` while preserving every segment in `segments`.
 * Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_cdt2d_build(const geomlib_vec2 *points, size_t count,
                        const uint32_t *segments, size_t segment_count,
                        geomlib_cdt2d *out);

/* Is (a,b) an edge of the triangulation? */
int geomlib_cdt2d_has_edge(const geomlib_cdt2d *c, uint32_t a, uint32_t b);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CDT2D_H */
