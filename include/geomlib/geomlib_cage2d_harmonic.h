/*
 * geomlib - 2D harmonic cage coordinates
 *
 * Computes positive harmonic-like coordinates inside a convex polygon using
 * inverse-distance weighted subtended angles, with affine precision.
 */

#ifndef GEOMLIB_CAGE2D_HARMONIC_H
#define GEOMLIB_CAGE2D_HARMONIC_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Weights are written to `weights` (n entries). Returns GEOMLIB_MESH_SUCCESS. */
int geomlib_cage2d_harmonic_weights(const geomlib_vec2 *cage, size_t n,
                                    geomlib_vec2 q, double *weights);

/* Deform query points from cage to deformed. */
int geomlib_cage2d_harmonic_deform(const geomlib_vec2 *cage,
                                   const geomlib_vec2 *deformed, size_t n,
                                   const geomlib_vec2 *points, size_t count,
                                   geomlib_vec2 *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CAGE2D_HARMONIC_H */
