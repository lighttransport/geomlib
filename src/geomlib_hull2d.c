/*
 * geomlib - 2D convex hull (see geomlib_hull2d.h)
 *
 * Clean-room monotone-chain construction. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_hull2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_hull2d geomlib_hull2d_create(void)
{
    geomlib_hull2d h;
    memset(&h, 0, sizeof(h));
    return h;
}

void geomlib_hull2d_free(geomlib_hull2d *h)
{
    if (!h) return;
    free(h->vertices);
    memset(h, 0, sizeof(*h));
}

static int cmp_point(const void *pa, const void *pb)
{
    const geomlib_vec2 *a=(const geomlib_vec2*)pa, *b=(const geomlib_vec2*)pb;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}

static double cross(const geomlib_vec2 *o, const geomlib_vec2 *a,
                    const geomlib_vec2 *b)
{
    return ((double)a->x-o->x)*((double)b->y-o->y)
         - ((double)a->y-o->y)*((double)b->x-o->x);
}

int geomlib_convex_hull2d(const geomlib_vec2 *points, size_t count,
                          geomlib_hull2d *out)
{
    geomlib_hull2d result = geomlib_hull2d_create();
    geomlib_vec2 *sorted=NULL, *hull=NULL;
    size_t n=0, k, i;

    if (!out || !points || !count) return GEOMLIB_MESH_ERROR_INVALID;
    sorted = (geomlib_vec2 *)malloc(count*sizeof(geomlib_vec2));
    hull = (geomlib_vec2 *)malloc((2*count+2)*sizeof(geomlib_vec2));
    if (!sorted || !hull) goto fail;

    for (i = 0; i < count; ++i) {
        size_t j;
        int dup = 0;
        for (j = 0; j < n; ++j)
            if (sorted[j].x==points[i].x && sorted[j].y==points[i].y) { dup=1; break; }
        if (!dup) sorted[n++] = points[i];
    }
    if (n < 1) goto fail;
    qsort(sorted, n, sizeof(geomlib_vec2), cmp_point);
    if (n <= 2) {
        result.vertices = (geomlib_vec2 *)malloc(n*sizeof(geomlib_vec2));
        if (!result.vertices) goto fail;
        memcpy(result.vertices, sorted, n*sizeof(geomlib_vec2));
        result.count = n;
        *out = result;
        free(sorted); free(hull);
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Lower chain. */
    k = 0;
    for (i = 0; i < n; ++i) {
        while (k >= 2 && cross(&hull[k-2], &hull[k-1], &sorted[i]) <= 0.0) --k;
        hull[k++] = sorted[i];
    }
    /* Upper chain. */
    {
        size_t lower = k+1;
        for (i = n; i-- > 0;) {
            while (k >= lower && cross(&hull[k-2], &hull[k-1], &sorted[i]) <= 0.0) --k;
            hull[k++] = sorted[i];
        }
        /* Last point repeats the first; drop it. */
        if (k > 1 && hull[k-1].x==hull[0].x && hull[k-1].y==hull[0].y) --k;
    }

    result.vertices = (geomlib_vec2 *)malloc((k?k:1)*sizeof(geomlib_vec2));
    if (!result.vertices) goto fail;
    memcpy(result.vertices, hull, k*sizeof(geomlib_vec2));
    result.count = k;
    *out = result;
    free(sorted); free(hull);
    return GEOMLIB_MESH_SUCCESS;

fail:
    free(sorted); free(hull);
    geomlib_hull2d_free(&result);
    return GEOMLIB_MESH_ERROR_MEMORY;
}
