/*
 * geomlib - Clip a triangle mesh by a plane and cap the cross-section
 *
 * Keeps one half-space of the mesh and closes the cut with a triangulated cap.
 */

#ifndef GEOMLIB_MESH_CLIP_H
#define GEOMLIB_MESH_CLIP_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Keep the side where sign(n·(p - point)) matches `keep_positive`. */
int geomlib_clip_mesh_plane(const geomlib_trimesh *mesh,
                            geomlib_vec3 plane_point, geomlib_vec3 plane_normal,
                            int keep_positive, geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MESH_CLIP_H */
