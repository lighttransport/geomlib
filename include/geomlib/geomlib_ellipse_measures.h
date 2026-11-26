/*
 * geomlib - Ellipse area and perimeter
 *
 * Area is exact; perimeter uses Ramanujan's second approximation.
 */

#ifndef GEOMLIB_ELLIPSE_MEASURES_H
#define GEOMLIB_ELLIPSE_MEASURES_H

#ifdef __cplusplus
extern "C" {
#endif

double geomlib_ellipse_area(float a, float b);
double geomlib_ellipse_perimeter(float a, float b);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ELLIPSE_MEASURES_H */
