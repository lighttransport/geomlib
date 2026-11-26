/*
 * geomlib - Scalar-field contour lines on a triangle mesh
 *
 * Extracts the zero/iso level set of a per-vertex field as surface segments.
 */

#ifndef GEOMLIB_MESH_CONTOUR_H
#define GEOMLIB_MESH_CONTOUR_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *segments;   /* two endpoints per segment */
    size_t segment_count;
} geomlib_mesh_contour;

geomlib_mesh_contour geomlib_mesh_contour_create(void);
void geomlib_mesh_contour_free(geomlib_mesh_contour *c);

/* Extract `iso` level of `values` (one per mesh vertex). */
int geomlib_mesh_contour_lines(const geomlib_trimesh *mesh,
                               const float *values, float iso,
                               geomlib_mesh_contour *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MESH_CONTOUR_H */
