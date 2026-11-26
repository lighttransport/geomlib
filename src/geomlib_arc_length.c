/*
 * geomlib - Polyline arc length (see geomlib_arc_length.h)
 *
 * Clean-room cumulative-distance table with linear interpolation. C11, no
 * external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_arc_length.h"

#include <math.h>

static double dist(const geomlib_vec3 *a, const geomlib_vec3 *b)
{
    double dx = (double)b->x - a->x, dy = (double)b->y - a->y,
           dz = (double)b->z - a->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

int geomlib_polyline_arc_table(const geomlib_vec3 *points, size_t n,
                               double *cumulative, double *total)
{
    size_t i;
    double acc = 0.0;
    if (!points || !cumulative || !n) return GEOMLIB_MESH_ERROR_INVALID;
    cumulative[0] = 0.0;
    for (i = 1; i < n; ++i) {
        acc += dist(&points[i - 1], &points[i]);
        cumulative[i] = acc;
    }
    if (total) *total = acc;
    return GEOMLIB_MESH_SUCCESS;
}

double geomlib_polyline_length(const geomlib_vec3 *points, size_t n)
{
    double total = 0.0;
    size_t i;
    if (!points || n < 2) return 0.0;
    for (i = 1; i < n; ++i) total += dist(&points[i - 1], &points[i]);
    return total;
}

geomlib_vec3 geomlib_polyline_at_length(const geomlib_vec3 *points, size_t n,
                                        const double *cumulative, double s)
{
    size_t i;
    if (!points || !cumulative || !n) return geomlib_vec3_create(0, 0, 0);
    if (n == 1) return points[0];
    if (s <= 0.0) return points[0];
    if (s >= cumulative[n - 1]) return points[n - 1];
    for (i = 1; i < n; ++i) {
        if (cumulative[i] >= s) {
            double seg = cumulative[i] - cumulative[i - 1];
            double t = seg > 1e-18 ? (s - cumulative[i - 1]) / seg : 0.0;
            return geomlib_vec3_add(
                points[i - 1],
                geomlib_vec3_mul(geomlib_vec3_sub(points[i], points[i - 1]),
                                 (float)t));
        }
    }
    return points[n - 1];
}
