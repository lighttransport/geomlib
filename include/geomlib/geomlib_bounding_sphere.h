/*
 * geomlib - Ritter bounding sphere
 *
 * Computes a fast enclosing sphere of a point set, with all points inside.
 */

#ifndef GEOMLIB_BOUNDING_SPHERE_H
#define GEOMLIB_BOUNDING_SPHERE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 center;
    float radius;
} geomlib_bounding_sphere;

/* Enclose all `count` points. Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_bounding_sphere_points(const geomlib_vec3 *points, size_t count,
                                   geomlib_bounding_sphere *out);

/* Enclose a triangle mesh's vertices. */
int geomlib_bounding_sphere_mesh(const geomlib_trimesh *mesh,
                                 geomlib_bounding_sphere *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_BOUNDING_SPHERE_H */
