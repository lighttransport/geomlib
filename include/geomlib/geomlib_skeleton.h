/*
 * geomlib - Curve skeleton via mean curvature flow
 *
 * Collapses a closed surface by mean-curvature flow with edge collapses that
 * preserve topology, leaving a 1D curve skeleton represented as vertices and
 * edges between contracted vertex clusters.
 */

#ifndef GEOMLIB_SKELETON_H
#define GEOMLIB_SKELETON_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *points;      /* Skeleton node positions */
    size_t point_count;
    uint32_t *edges;           /* 2 indices per edge */
    size_t edge_count;
    uint32_t *original_vertices; /* Representative original vertex per node */
} geomlib_curve_skeleton;

geomlib_curve_skeleton geomlib_curve_skeleton_create(void);
void geomlib_curve_skeleton_free(geomlib_curve_skeleton *s);

typedef struct {
    float collapse_ratio;   /* Stop when vertex count drops below ratio (default 0.02) */
    int max_iterations;     /* Default 200 */
    float feature_angle;    /* Protect sharp edges (degrees, default 45) */
} geomlib_skeleton_params;

geomlib_skeleton_params geomlib_skeleton_params_default(void);

/*
 * Extract a curve skeleton from a closed triangle mesh. The input is not
 * modified. Returns the contracted graph.
 */
int geomlib_extract_curve_skeleton(const geomlib_trimesh *mesh,
                                   const geomlib_skeleton_params *params,
                                   geomlib_curve_skeleton *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SKELETON_H */
