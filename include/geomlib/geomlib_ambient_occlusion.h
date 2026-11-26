/*
 * geomlib - Vertex ambient occlusion by hemisphere ray casting
 *
 * Estimates the fraction of the outward hemisphere visible from each vertex.
 */

#ifndef GEOMLIB_AMBIENT_OCCLUSION_H
#define GEOMLIB_AMBIENT_OCCLUSION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Write one occlusion value in [0,1] per vertex (1 = fully open). */
int geomlib_mesh_ambient_occlusion(geomlib_trimesh *mesh, int samples,
                                   float max_distance, float *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_AMBIENT_OCCLUSION_H */
