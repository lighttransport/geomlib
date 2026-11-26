/*
 * geomlib - Ellipse measures (see geomlib_ellipse_measures.h)
 *
 * Clean-room area and Ramanujan perimeter. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_ellipse_measures.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double geomlib_ellipse_area(float a, float b)
{
    return M_PI * fabs((double)a) * fabs((double)b);
}

double geomlib_ellipse_perimeter(float a, float b)
{
    double aa = fabs((double)a), bb = fabs((double)b);
    double sum = aa + bb;
    double h, d;
    if (!(sum > 1e-12)) return 0.0;
    h = (aa - bb) * (aa - bb) / (sum * sum);
    d = sqrt(fmax(0.0, 4.0 - 3.0 * h));
    return M_PI * sum * (1.0 + 3.0 * h / (10.0 + d));
}
