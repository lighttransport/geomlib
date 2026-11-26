/*
 * geomlib - 2D polygon medial axis
 *
 * Approximates the medial axis with Voronoi vertices of densely sampled
 * boundary points. A Voronoi vertex is kept when it lies strictly inside the
 * polygon and its nearest samples belong to at least two different edges.
 */

#ifndef GEOMLIB_MEDIAL2D_H
#define GEOMLIB_MEDIAL2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *segments;   /* two endpoints per segment */
    size_t segment_count;
} geomlib_medial2d;

geomlib_medial2d geomlib_medial2d_create(void);
void geomlib_medial2d_free(geomlib_medial2d *m);

/* Build the medial axis of a simple polygon. `spacing` <= 0 selects
 * 0.02 * bounding-box diagonal. */
int geomlib_medial_axis2d(const geomlib_vec2 *polygon, size_t count,
                          float spacing, geomlib_medial2d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MEDIAL2D_H */
