/*
 * geomlib - Cross fields on triangle meshes
 *
 * Generates a 4-directional (cross) field from principal curvature directions
 * and combs it across the surface so neighbouring directions differ minimally
 * modulo 90 degrees.
 */

#ifndef GEOMLIB_CROSS_FIELD_H
#define GEOMLIB_CROSS_FIELD_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float *angles;      /* Cross-field direction per vertex, in local frame */
    float *normals_x;   /* Local frame first tangent per vertex */
    float *normals_y;   /* Local frame second tangent per vertex */
    size_t count;
} geomlib_cross_field;

geomlib_cross_field geomlib_cross_field_create(size_t count);
void geomlib_cross_field_free(geomlib_cross_field *f);

/* Generate and comb a cross field. */
int geomlib_compute_cross_field(const geomlib_trimesh *mesh,
                                geomlib_cross_field *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CROSS_FIELD_H */
