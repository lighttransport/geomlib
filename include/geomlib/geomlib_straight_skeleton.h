/*
 * geomlib - 2D straight skeleton of a simple polygon
 *
 * Computes the straight skeleton via the "wavefront propagation" approach:
 * at each event, adjacent edges are offset inward at unit speed and their
 * intersections are advanced; split and vertex events are handled by
 * collapsing zero-length edges. The output is the set of skeleton segments
 * with their associated collapse times.
 */

#ifndef GEOMLIB_STRAIGHT_SKELETON_H
#define GEOMLIB_STRAIGHT_SKELETON_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *points;      /* Skeleton nodes (x, y) */
    size_t point_count;
    uint32_t *edges;           /* 2 indices per skeleton edge */
    size_t edge_count;
    float *times;              /* Collapse time per skeleton node */
} geomlib_straight_skeleton;

geomlib_straight_skeleton geomlib_straight_skeleton_create(void);
void geomlib_straight_skeleton_free(geomlib_straight_skeleton *s);

/*
 * Compute the interior straight skeleton of a simple polygon in 2D.
 * Points are given in CCW order; the polygon may be non-convex.
 */
int geomlib_compute_straight_skeleton(const geomlib_polygon2d *polygon,
                                      geomlib_straight_skeleton *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_STRAIGHT_SKELETON_H */
