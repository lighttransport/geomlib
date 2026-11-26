/*
 * geomlib - Algebraic least-squares 2D circle fitting
 *
 * Fits a circle to 2D points using a linear least-squares solve on the
 * circle equation.
 */

#ifndef GEOMLIB_CIRCLE_FIT_H
#define GEOMLIB_CIRCLE_FIT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 center;
    float radius;
    float rms;
} geomlib_circle_fit;

int geomlib_fit_circle(const geomlib_vec2 *points, size_t count,
                       geomlib_circle_fit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_FIT_H */
