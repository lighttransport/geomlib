/*
 * geomlib - Direct least-squares ellipse fitting
 *
 * Fits a general conic with the ellipse constraint to 2D points using
 * Fitzgibbon's scatter-matrix eigen formulation.
 */

#ifndef GEOMLIB_ELLIPSE_FIT_H
#define GEOMLIB_ELLIPSE_FIT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 center;
    float axis_a, axis_b;   /* semi-axes (a >= b) */
    float angle;            /* major-axis angle in radians */
    float rms;
} geomlib_ellipse_fit;

int geomlib_fit_ellipse(const geomlib_vec2 *points, size_t count,
                        geomlib_ellipse_fit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ELLIPSE_FIT_H */
