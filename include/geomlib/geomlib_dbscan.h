/*
 * geomlib - DBSCAN density-based point-cloud clustering
 *
 * Labels dense neighbourhoods as clusters and marks sparse points as noise.
 */

#ifndef GEOMLIB_DBSCAN_H
#define GEOMLIB_DBSCAN_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GEOMLIB_DBSCAN_NOISE UINT32_MAX

typedef struct {
    float eps;            /* neighbourhood radius (0 = 2 * mean spacing) */
    int min_points;       /* minimum neighbours including the point (default 5) */
} geomlib_dbscan_params;

geomlib_dbscan_params geomlib_dbscan_params_default(void);

/* Write one label per point; noise points get GEOMLIB_DBSCAN_NOISE. */
int geomlib_pointcloud_dbscan(const geomlib_pointcloud *cloud,
                              const geomlib_dbscan_params *params,
                              uint32_t *labels, size_t *cluster_count);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DBSCAN_H */
