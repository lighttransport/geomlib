/*
 * geomlib - Minimum distance between two triangle meshes
 *
 * Finds the closest pair of points by querying each mesh's vertices against
 * the other mesh's BVH, then refining from the better seed.
 */

#ifndef GEOMLIB_MESH_DISTANCE_H
#define GEOMLIB_MESH_DISTANCE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 point_a, point_b;
    float distance;
} geomlib_mesh_distance;

/* Minimum Euclidean distance between the surfaces of `a` and `b`. */
int geomlib_mesh_mesh_distance(const geomlib_trimesh *a,
                               const geomlib_trimesh *b,
                               geomlib_mesh_distance *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MESH_DISTANCE_H */
