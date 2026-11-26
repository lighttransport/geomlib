/*
 * geomlib - Circle areas (see geomlib_circle_area.h)
 *
 * Clean-room analytic lens formula. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_circle_area.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double geomlib_circle_intersection_area(double r0, double r1, double d)
{
    double a = fabs(r0), b = fabs(r1), dist = fabs(d);
    double a2, b2, d2, c0, c1, term;

    if (!(a > 0.0 && b > 0.0)) return 0.0;
    if (dist >= a + b) return 0.0;
    if (dist <= fabs(a - b)) {
        double m = a < b ? a : b;
        return M_PI * m * m;
    }
    a2 = a * a; b2 = b * b; d2 = dist * dist;
    c0 = (d2 + a2 - b2) / (2.0 * dist * a);
    c1 = (d2 + b2 - a2) / (2.0 * dist * b);
    if (c0 > 1.0) c0 = 1.0;
    if (c0 < -1.0) c0 = -1.0;
    if (c1 > 1.0) c1 = 1.0;
    if (c1 < -1.0) c1 = -1.0;
    term = (-dist + a + b) * (dist + a - b) * (dist - a + b) * (dist + a + b);
    return a2 * acos(c0) + b2 * acos(c1) - 0.5 * sqrt(fmax(0.0, term));
}

double geomlib_circle_union_area(double r0, double r1, double d)
{
    return M_PI * r0 * r0 + M_PI * r1 * r1
         - geomlib_circle_intersection_area(r0, r1, d);
}
