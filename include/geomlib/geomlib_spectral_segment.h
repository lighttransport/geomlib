/*
 * geomlib - Spectral mesh segmentation
 *
 * Embeds vertices with the lowest nontrivial Laplace-Beltrami eigenvectors
 * and partitions them with deterministic k-means.
 */

#ifndef GEOMLIB_SPECTRAL_SEGMENT_H
#define GEOMLIB_SPECTRAL_SEGMENT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Label each vertex into `k` clusters. `eigen_count` should be >= k. */
int geomlib_spectral_segment(const geomlib_trimesh *mesh, size_t k,
                             size_t eigen_count, uint32_t *labels);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPECTRAL_SEGMENT_H */
