/*
 * geomlib - Visvalingam-Whyatt simplification (see geomlib_vw.h)
 *
 * Clean-room repeated-minimum effective-area removal. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_vw.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

static double tri_area2(const geomlib_vec3 *points, size_t i, size_t j, size_t k)
{
    double ax = (double)points[j].x - points[i].x;
    double ay = (double)points[j].y - points[i].y;
    double az = (double)points[j].z - points[i].z;
    double bx = (double)points[k].x - points[i].x;
    double by = (double)points[k].y - points[i].y;
    double bz = (double)points[k].z - points[i].z;
    double cx = ay * bz - az * by;
    double cy = az * bx - ax * bz;
    double cz = ax * by - ay * bx;
    return sqrt(cx * cx + cy * cy + cz * cz);
}

double geomlib_vw_effective_area(const geomlib_vec3 *points, size_t count,
                                 size_t i)
{
    if (!points || i == 0 || i + 1 >= count) return 0.0;
    return 0.5 * tri_area2(points, i - 1, i, i + 1);
}

int geomlib_simplify_polyline_vw(const geomlib_vec3 *points, size_t count,
                                 double area_epsilon, uint32_t *indices,
                                 size_t *out_count)
{
    size_t *prev = NULL, *next = NULL, *alive = NULL;
    size_t i, n;

    if (!points || !indices || !out_count || count < 2 || !(area_epsilon >= 0.0))
        return GEOMLIB_MESH_ERROR_INVALID;
    prev = (size_t *)malloc(count * sizeof(size_t));
    next = (size_t *)malloc(count * sizeof(size_t));
    alive = (size_t *)malloc(count * sizeof(size_t));
    if (!prev || !next || !alive) goto fail;
    for (i = 0; i < count; ++i) {
        prev[i] = i ? i - 1 : i;
        next[i] = (i + 1 < count) ? i + 1 : i;
        alive[i] = 1;
    }
    n = count;

    for (;;) {
        size_t best = (size_t)-1;
        double best_area = DBL_MAX;
        for (i = 1; i + 1 < count; ++i) {
            double a;
            if (!alive[i]) continue;
            a = 0.5 * tri_area2(points, prev[i], i, next[i]);
            if (a < best_area) { best_area = a; best = i; }
        }
        if (best == (size_t)-1 || best_area > area_epsilon) break;

        alive[best] = 0;
        next[prev[best]] = next[best];
        prev[next[best]] = prev[best];
        n--;
        if (n <= 2) break;
    }

    {
        size_t k = 0;
        i = 0;
        while (k < n) {
            indices[k++] = (uint32_t)i;
            i = next[i];
        }
    }
    *out_count = n;
    free(prev); free(next); free(alive);
    return GEOMLIB_MESH_SUCCESS;

fail:
    free(prev); free(next); free(alive);
    return GEOMLIB_MESH_ERROR_MEMORY;
}
