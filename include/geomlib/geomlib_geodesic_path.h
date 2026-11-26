/*
 * geomlib - Shortest-path polyline between mesh vertices
 *
 * Reconstructs the edge-graph shortest path from source to target.
 */

#ifndef GEOMLIB_GEODESIC_PATH_H
#define GEOMLIB_GEODESIC_PATH_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *vertices;   /* path endpoints included, in traversal order */
    size_t count;
    float length;
} geomlib_geodesic_path;

geomlib_geodesic_path geomlib_geodesic_path_create(void);
void geomlib_geodesic_path_free(geomlib_geodesic_path *p);

int geomlib_geodesic_path_between(const geomlib_trimesh *mesh,
                                  uint32_t source, uint32_t target,
                                  geomlib_geodesic_path *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_PATH_H */
