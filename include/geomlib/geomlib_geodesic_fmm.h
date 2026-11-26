/*
 * geomlib - Fast-marching-style surface geodesics
 *
 * Computes approximate geodesic distances on a triangle mesh from a source
 * vertex: an edge-graph Dijkstra initialization followed by iterative
 * triangle-interior relaxation (unfolded wavefront propagation).
 */

#ifndef GEOMLIB_GEODESIC_FMM_H
#define GEOMLIB_GEODESIC_FMM_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float *distances;
    size_t count;
} geomlib_geodesic_fmm_result;

geomlib_geodesic_fmm_result geomlib_geodesic_fmm_result_create(size_t count);
void geomlib_geodesic_fmm_result_free(geomlib_geodesic_fmm_result *r);

/* Geodesic distances from `source` on a closed triangle mesh. */
int geomlib_geodesic_fast_marching(const geomlib_trimesh *mesh,
                                   uint32_t source,
                                   int relaxation_passes,
                                   geomlib_geodesic_fmm_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_FMM_H */
