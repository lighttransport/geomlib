/*
 * geomlib - Mesh saliency and saliency-based segmentation
 *
 * Saliency measures local "interestingness" of a vertex from smoothed
 * curvature at multiple spatial scales (Lee et al. 2005). High-saliency
 * regions are then segmented by flood fill for shape analysis.
 */

#ifndef GEOMLIB_SALIENCY_H
#define GEOMLIB_SALIENCY_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float *values;        /* Normalized saliency per vertex, length vertex_count */
    size_t vertex_count;
} geomlib_saliency_result;

geomlib_saliency_result geomlib_saliency_result_create(size_t vertex_count);
void geomlib_saliency_result_free(geomlib_saliency_result *r);

/*
 * Compute multi-scale curvature saliency. `scales` is the number of
 * smoothing scales (default 3 when 0).
 */
int geomlib_compute_saliency(const geomlib_trimesh *mesh,
                             int scales,
                             geomlib_saliency_result *out);

/* Segment vertices into `region_count` regions seeded at the most salient
 * vertices. Returns region labels per vertex (0..region_count-1). */
int geomlib_segment_by_saliency(const geomlib_trimesh *mesh,
                                const geomlib_saliency_result *saliency,
                                size_t region_count,
                                uint32_t *labels);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SALIENCY_H */
