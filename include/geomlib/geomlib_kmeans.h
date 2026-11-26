/*
 * geomlib - Deterministic k-means clustering for point clouds
 *
 * Farthest-point seeding followed by Lloyd iterations.
 */

#ifndef GEOMLIB_KMEANS_H
#define GEOMLIB_KMEANS_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Assign each point to one of k clusters. `labels` has cloud->point_count
 * entries; `centroids` has k entries (may be NULL). Returns cluster count. */
int geomlib_pointcloud_kmeans(const geomlib_pointcloud *cloud, size_t k,
                              int max_iterations,
                              uint32_t *labels, geomlib_vec3 *centroids);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_KMEANS_H */
