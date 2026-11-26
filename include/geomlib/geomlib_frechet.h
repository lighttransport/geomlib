/*
 * geomlib - Discrete Fréchet distance between polylines
 *
 * Computes the classic discrete Fréchet distance via dynamic programming.
 */

#ifndef GEOMLIB_FRECHET_H
#define GEOMLIB_FRECHET_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Discrete Fréchet distance between polylines a (na points) and b (nb). */
double geomlib_frechet_distance(const geomlib_vec3 *a, size_t na,
                                const geomlib_vec3 *b, size_t nb);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_FRECHET_H */
