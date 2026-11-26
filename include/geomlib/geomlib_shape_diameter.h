/*
 * geomlib - Shape Diameter Function
 *
 * For each surface vertex, cast a cone of rays around the outward normal and
 * summarize the ray lengths into a scale-invariant local thickness measure.
 */

#ifndef GEOMLIB_SHAPE_DIAMETER_H
#define GEOMLIB_SHAPE_DIAMETER_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float *values;       /* Normalized shape diameter per vertex */
    size_t count;
} geomlib_shape_diameter;

geomlib_shape_diameter geomlib_shape_diameter_create(size_t count);
void geomlib_shape_diameter_free(geomlib_shape_diameter *d);

/* Compute the shape diameter function of a closed mesh. */
int geomlib_compute_shape_diameter(const geomlib_trimesh *mesh,
                                   int ray_count,
                                   geomlib_shape_diameter *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SHAPE_DIAMETER_H */
