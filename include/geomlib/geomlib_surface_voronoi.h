/*
 * geomlib - Surface Voronoi segmentation
 *
 * Assigns every mesh vertex to its nearest site in geodesic (mesh-edge)
 * distance, producing a partition of the surface into Voronoi regions.
 */

#ifndef GEOMLIB_SURFACE_VORONOI_H
#define GEOMLIB_SURFACE_VORONOI_H

#include "geomlib_geodesic.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *labels;      /* one nearest-site index per vertex */
    size_t *sizes;         /* region sizes, one per site */
    size_t vertex_count;
    size_t site_count;
} geomlib_surface_voronoi;

geomlib_surface_voronoi geomlib_surface_voronoi_create(void);
void geomlib_surface_voronoi_free(geomlib_surface_voronoi *v);

/* Partition `mesh` into `nsites` geodesic Voronoi regions. */
int geomlib_surface_voronoi_compute(const geomlib_trimesh *mesh,
                                    const uint32_t *sites, size_t nsites,
                                    geomlib_surface_voronoi *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SURFACE_VORONOI_H */
