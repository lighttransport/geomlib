/*
 * geomlib - Circle intersection and union areas
 *
 * Analytic lens area of two circles in the plane.
 */

#ifndef GEOMLIB_CIRCLE_AREA_H
#define GEOMLIB_CIRCLE_AREA_H

#ifdef __cplusplus
extern "C" {
#endif

/* Area common to two circles of radii r0, r1 whose centres are distance d
 * apart. Returns 0 for disjoint/tangent circles and pi*min(r0,r1)^2 for
 * containment. */
double geomlib_circle_intersection_area(double r0, double r1, double d);

/* Area covered by either circle. */
double geomlib_circle_union_area(double r0, double r1, double d);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_AREA_H */
