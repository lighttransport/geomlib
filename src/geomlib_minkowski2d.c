/*
 * geomlib - 2D Minkowski sum (see geomlib_minkowski2d.h)
 *
 * For convex polygons, Minkowski(A,B) is the convex hull of the na*nb
 * pairwise vertex sums. This clean-room construction reuses the monotone-chain
 * hull and is exact for convex inputs.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_minkowski2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_minkowski_sum2d(const geomlib_vec2 *a, size_t na,
                            const geomlib_vec2 *b, size_t nb,
                            geomlib_hull2d *out)
{
    geomlib_vec2 *sums;
    geomlib_hull2d result;
    size_t i, j, k = 0;

    if (!out || !a || !b || !na || !nb) return GEOMLIB_MESH_ERROR_INVALID;
    if (na > SIZE_MAX/nb) return GEOMLIB_MESH_ERROR_MEMORY;
    sums = (geomlib_vec2 *)malloc(na*nb*sizeof(geomlib_vec2));
    if (!sums) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i = 0; i < na; ++i)
        for (j = 0; j < nb; ++j)
            sums[k++] = geomlib_vec2_create(a[i].x+b[j].x, a[i].y+b[j].y);

    result = geomlib_hull2d_create();
    if (geomlib_convex_hull2d(sums, k, &result) != GEOMLIB_MESH_SUCCESS) {
        free(sums);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    free(sums);
    *out = result;
    return GEOMLIB_MESH_SUCCESS;
}
