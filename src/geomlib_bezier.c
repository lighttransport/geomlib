/*
 * geomlib - Bezier curves (see geomlib_bezier.h)
 *
 * Clean-room de Casteljau algorithm. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_bezier.h"

#include <math.h>
#include <stdlib.h>

static double clamp01(double t)
{
    if (t < 0.0) return 0.0;
    if (t > 1.0) return 1.0;
    return t;
}

/* Build the full de Casteljau triangle: level[j] holds n-j points. */
static int casteljau(const geomlib_vec3 *control, size_t n, double t,
                     double (*levels)[3])
{
    size_t i, j, k;
    for (i = 0; i < n; ++i) {
        levels[i][0] = (double)control[i].x;
        levels[i][1] = (double)control[i].y;
        levels[i][2] = (double)control[i].z;
    }
    for (j = 1; j < n; ++j)
        for (i = 0; i + j < n; ++i) {
            double *a = levels[(j - 1) * n + i];
            double *b = levels[(j - 1) * n + i + 1];
            double *d = levels[j * n + i];
            for (k = 0; k < 3; ++k) d[k] = (1.0 - t) * a[k] + t * b[k];
        }
    return 1;
}

static geomlib_vec3 from_double(const double *v)
{
    return geomlib_vec3_create((float)v[0], (float)v[1], (float)v[2]);
}

geomlib_vec3 geomlib_bezier_evaluate(const geomlib_vec3 *control, size_t n,
                                     double t)
{
    double (*levels)[3];
    geomlib_vec3 r;
    if (!control || n < 1) return geomlib_vec3_create(0, 0, 0);
    if (n == 1) return control[0];
    levels = (double (*)[3])malloc(n * n * 3 * sizeof(double));
    if (!levels) return control[0];
    casteljau(control, n, clamp01(t), levels);
    r = from_double(levels[(n - 1) * n]);
    free(levels);
    return r;
}

geomlib_vec3 geomlib_bezier_derivative(const geomlib_vec3 *control, size_t n,
                                       double t)
{
    geomlib_vec3 *diff;
    geomlib_vec3 r;
    size_t i;
    if (!control || n < 2) return geomlib_vec3_create(0, 0, 0);
    diff = (geomlib_vec3 *)malloc((n - 1) * sizeof(geomlib_vec3));
    if (!diff) return geomlib_vec3_create(0, 0, 0);
    for (i = 0; i + 1 < n; ++i)
        diff[i] = geomlib_vec3_mul(
            geomlib_vec3_sub(control[i + 1], control[i]), (float)(n - 1));
    r = geomlib_bezier_evaluate(diff, n - 1, t);
    free(diff);
    return r;
}

void geomlib_bezier_subdivide(const geomlib_vec3 *control, size_t n, double t,
                              geomlib_vec3 *left, geomlib_vec3 *right)
{
    double (*levels)[3];
    size_t j;
    if (!control || !left || !right || n < 1) return;
    if (n == 1) { left[0] = right[0] = control[0]; return; }
    levels = (double (*)[3])malloc(n * n * 3 * sizeof(double));
    if (!levels) return;
    casteljau(control, n, clamp01(t), levels);
    for (j = 0; j < n; ++j) {
        left[j] = from_double(levels[j * n]);               /* Q[j][0] */
        right[j] = from_double(levels[(n - 1 - j) * n + j]); /* Q[n-1-j][j] */
    }
    free(levels);
}
