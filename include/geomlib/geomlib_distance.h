/*
 * geomlib - Mesh distance queries with error bounds
 */

#ifndef GEOMLIB_DISTANCE_H
#define GEOMLIB_DISTANCE_H

#include "geomlib_mesh.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float lower_bound;  /* largest sampled point-to-surface distance */
    float upper_bound;  /* rigorous upper bound: lower + sampling error */
    size_t samples;     /* number of sampled surface points per direction */
} geomlib_hausdorff_result;

/*
 * Bounded-error symmetric Hausdorff distance between two triangle meshes.
 *
 * Because the point-to-surface distance function is 1-Lipschitz, a point in a
 * source triangle is at distance at most (minimum vertex distance + longest
 * edge length) from the target surface. Vertices, edge midpoints and centroids
 * are sampled for a tight lower bound, and that inequality supplies a rigorous
 * upper bound. Both directions are evaluated.
 */
int geomlib_hausdorff_bounded(const geomlib_trimesh *a,
                              const geomlib_trimesh *b,
                              geomlib_hausdorff_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DISTANCE_H */
