/*
 * geomlib - Farthest-point sampling on surface vertices
 *
 * Greedily selects the vertex with the largest distance to the current sample
 * set, producing well-spread samples with a monotonically decreasing covering
 * radius. Distances are geodesic, with a discrete edge-propagation fallback.
 */

#ifndef GEOMLIB_FPS_H
#define GEOMLIB_FPS_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *indices;         /* k sampled vertex indices */
    float *radii;              /* covering radius after each added sample */
    size_t count;
} geomlib_fps_result;

geomlib_fps_result geomlib_fps_result_create(void);
void geomlib_fps_result_free(geomlib_fps_result *r);

/* Sample k vertices, starting from `seed`. Returns GEOMLIB_MESH_SUCCESS. */
int geomlib_farthest_point_sampling(const geomlib_trimesh *mesh, size_t k,
                                    uint32_t seed, geomlib_fps_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_FPS_H */
