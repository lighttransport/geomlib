/*
 * geomlib - Algebraic least-squares 3D sphere fitting
 *
 * Fits a sphere to 3D points using a linear least-squares solve on the
 * sphere equation.
 */

#ifndef GEOMLIB_SPHERE_FIT_H
#define GEOMLIB_SPHERE_FIT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 center;
    float radius;
    float rms;
} geomlib_sphere_fit;

int geomlib_fit_sphere(const geomlib_vec3 *points, size_t count,
                       geomlib_sphere_fit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPHERE_FIT_H */
