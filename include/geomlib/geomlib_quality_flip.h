/*
 * geomlib - Surface mesh triangle-quality improvement by edge flips
 *
 * Flips interior edges when the resulting quadrilateral is convex and the
 * minimum triangle angle improves. Boundary edges are never flipped.
 */

#ifndef GEOMLIB_QUALITY_FLIP_H
#define GEOMLIB_QUALITY_FLIP_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Improve triangle quality in place. Returns the number of flips performed. */
int geomlib_improve_mesh_quality(geomlib_trimesh *mesh, int max_passes);

/* Minimum triangle angle over the mesh, in radians. */
float geomlib_mesh_min_angle(const geomlib_trimesh *mesh);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_QUALITY_FLIP_H */
