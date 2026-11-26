/*
 * geomlib - Robust geometric predicates
 *
 * Clean-room implementation of adaptive exact predicates based on the
 * published expansion-arithmetic method.
 *
 * Reference (paper only; no target library source was consulted):
 *   Jonathan R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and
 *   Fast Robust Geometric Predicates", Discrete & Computational Geometry,
 *   vol. 18, no. 3, pp. 305-363, 1997.
 *
 * All predicates return an exact sign for finite IEEE-754 double inputs:
 *   > 0  positive orientation / point inside circle or sphere
 *   = 0  exactly degenerate
 *   < 0  negative orientation / point outside circle or sphere
 *
 * A fast floating-point determinant with a conservative error bound is used
 * as a filter; only when the sign is uncertain is an exact expansion computed.
 * No long double, no compiler extensions: safe for C11 and WebAssembly.
 */

#ifndef GEOMLIB_PREDICATES_H
#define GEOMLIB_PREDICATES_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sign of the 2x2 orientation determinant of points a, b, c:
 *   (b - a) x (c - a)
 * Positive for counter-clockwise order. */
int geomlib_pred_orient2d(double ax, double ay,
                          double bx, double by,
                          double cx, double cy);

/* Sign of the 3x3 orientation determinant of points a, b, c (2D float). */
int geomlib_pred_orient2df(float ax, float ay,
                           float bx, float by,
                           float cx, float cy);

/* Sign of the 3x3 orientation determinant of points a, b, c in 3D:
 *   det [b-a, c-a, d-a]  (points a, b, c, d) */
int geomlib_pred_orient3d(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double cx, double cy, double cz,
                          double dx, double dy, double dz);

/* Sign of the in-circle determinant for points a, b, c (assumed CCW) and d:
 *   > 0 if d lies inside the circumcircle of a, b, c
 *   = 0 if cocircular
 *   < 0 if outside */
int geomlib_pred_incircle(double ax, double ay,
                          double bx, double by,
                          double cx, double cy,
                          double dx, double dy);

/* Sign of the in-sphere determinant for points a, b, c, d (assumed positively
 * oriented) and query e:
 *   > 0 if e lies inside the circumsphere of a, b, c, d
 *   = 0 if cospherical
 *   < 0 if outside */
int geomlib_pred_insphere(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double cx, double cy, double cz,
                          double dx, double dy, double dz,
                          double ex, double ey, double ez);

/* Generic exact sign of the determinant of a row-major n x n double matrix.
 * n must be in [1, 5]. Used internally by the predicates above and exposed
 * for callers that need an exact small determinant. */
int geomlib_pred_det_sign(int n, const double *matrix);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PREDICATES_H */
