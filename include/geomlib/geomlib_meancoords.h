/*
 * geomlib - 2D mean-value coordinate cage deformation
 *
 * Deforms points inside a polygon cage using positive, affine-precise
 * mean-value weights derived from the angles subtended by cage edges.
 */

#ifndef GEOMLIB_MEANCOORDS_H
#define GEOMLIB_MEANCOORDS_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Deform `count` query points. `cage` and `deformed` each hold `n` vertices.
 * Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_deform_mean_value2d(const geomlib_vec2 *cage,
                                const geomlib_vec2 *deformed, size_t n,
                                const geomlib_vec2 *points, size_t count,
                                geomlib_vec2 *out);

/* Mean-value weights of a single point (out array of n doubles). */
int geomlib_mean_value_weights2d(const geomlib_vec2 *cage, size_t n,
                                 geomlib_vec2 q, double *weights);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MEANCOORDS_H */
