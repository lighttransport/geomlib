/*
 * geomlib - Bezier curve evaluation and subdivision
 *
 * De Casteljau evaluation, subdivision and derivatives of arbitrary-degree
 * Bezier curves in 3D.
 */

#ifndef GEOMLIB_BEZIER_H
#define GEOMLIB_BEZIER_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Evaluate the Bezier curve with `n` control points at t (clamped to [0,1]). */
geomlib_vec3 geomlib_bezier_evaluate(const geomlib_vec3 *control, size_t n,
                                     double t);

/* First derivative at t. */
geomlib_vec3 geomlib_bezier_derivative(const geomlib_vec3 *control, size_t n,
                                       double t);

/* Split at t into two curves, each with `n` control points. The last point of
 * `left` and the first point of `right` are both the split point. */
void geomlib_bezier_subdivide(const geomlib_vec3 *control, size_t n, double t,
                              geomlib_vec3 *left, geomlib_vec3 *right);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_BEZIER_H */
