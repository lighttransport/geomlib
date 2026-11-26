/*
 * geomlib - Spatial sorting (see geomlib_spatial_sort.h)
 *
 * Clean-room bit-interleaving for Morton codes and the standard recursive
 * rotation Hilbert curve. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_spatial_sort.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static uint64_t part1by1(uint64_t x)
{
    x &= 0xffffffffULL;
    x = (x | (x<<16)) & 0x0000ffff0000ffffULL;
    x = (x | (x<<8))  & 0x00ff00ff00ff00ffULL;
    x = (x | (x<<4))  & 0x0f0f0f0f0f0f0f0fULL;
    x = (x | (x<<2))  & 0x3333333333333333ULL;
    x = (x | (x<<1))  & 0x5555555555555555ULL;
    return x;
}

static uint64_t compact1by1(uint64_t x)
{
    x &= 0x5555555555555555ULL;
    x = (x | (x>>1))  & 0x3333333333333333ULL;
    x = (x | (x>>2))  & 0x0f0f0f0f0f0f0f0fULL;
    x = (x | (x>>4))  & 0x00ff00ff00ff00ffULL;
    x = (x | (x>>8))  & 0x0000ffff0000ffffULL;
    x = (x | (x>>16)) & 0x00000000ffffffffULL;
    return x;
}

static uint64_t part1by2(uint64_t x)
{
    x &= 0xffffffffULL;
    x = (x | (x<<32)) & 0x1f00000000ffffffULL;
    x = (x | (x<<16)) & 0x1f0000ff0000ffULL;
    x = (x | (x<<8))  & 0x100f00f00f00f00fULL;
    x = (x | (x<<4))  & 0x10c30c30c30c30c3ULL;
    x = (x | (x<<2))  & 0x1249249249249249ULL;
    return x;
}

static uint64_t compact1by2(uint64_t x)
{
    x &= 0x1249249249249249ULL;
    x = (x | (x>>2))  & 0x10c30c30c30c30c3ULL;
    x = (x | (x>>4))  & 0x100f00f00f00f00fULL;
    x = (x | (x>>8))  & 0x1f0000ff0000ffULL;
    x = (x | (x>>16)) & 0x1f00000000ffffffULL;
    x = (x | (x>>32)) & 0x00000000ffffffffULL;
    return x;
}

uint64_t geomlib_morton2_encode(uint32_t x, uint32_t y)
{
    return part1by1(x) | (part1by1(y)<<1);
}

uint64_t geomlib_morton3_encode(uint32_t x, uint32_t y, uint32_t z)
{
    return part1by2(x) | (part1by2(y)<<1) | (part1by2(z)<<2);
}

void geomlib_morton2_decode(uint64_t code, uint32_t *x, uint32_t *y)
{
    if (x) *x = (uint32_t)compact1by1(code);
    if (y) *y = (uint32_t)compact1by1(code>>1);
}

void geomlib_morton3_decode(uint64_t code, uint32_t *x, uint32_t *y, uint32_t *z)
{
    if (x) *x = (uint32_t)compact1by2(code);
    if (y) *y = (uint32_t)compact1by2(code>>1);
    if (z) *z = (uint32_t)compact1by2(code>>2);
}

static void hilbert_rotate(uint32_t s, uint32_t *x, uint32_t *y,
                           uint32_t rx, uint32_t ry)
{
    if (ry == 0) {
        if (rx == 1) { *x = s-1-*x; *y = s-1-*y; }
        { uint32_t t = *x; *x = *y; *y = t; }
    }
}

uint64_t geomlib_hilbert2_encode(uint32_t side, uint32_t x, uint32_t y)
{
    uint64_t d = 0;
    uint32_t s;
    if (side < 2) return 0;
    for (s = side/2; s > 0; s /= 2) {
        uint32_t rx = (x & s) ? 1u : 0u;
        uint32_t ry = (y & s) ? 1u : 0u;
        d += (uint64_t)s*s*((3u*rx)^ry);
        hilbert_rotate(s, &x, &y, rx, ry);
    }
    return d;
}

void geomlib_hilbert2_decode(uint32_t side, uint64_t d,
                             uint32_t *outx, uint32_t *outy)
{
    uint32_t x = 0, y = 0, s;
    if (side < 2) { if (outx) *outx = 0; if (outy) *outy = 0; return; }
    for (s = 1; s < side; s *= 2) {
        uint32_t rx = (uint32_t)(1u & (d/2));
        uint32_t ry = (uint32_t)(1u & (d^rx));
        hilbert_rotate(s, &x, &y, rx, ry);
        x += s*rx;
        y += s*ry;
        d /= 4;
    }
    if (outx) *outx = x;
    if (outy) *outy = y;
}

typedef struct { uint64_t code; uint32_t id; } keyed;

static int cmp_keyed(const void *pa, const void *pb)
{
    const keyed *a=(const keyed*)pa, *b=(const keyed*)pb;
    if (a->code < b->code) return -1;
    if (a->code > b->code) return 1;
    return a->id < b->id ? -1 : (a->id > b->id ? 1 : 0);
}

static uint32_t quantize(double v, double lo, double hi, uint32_t levels)
{
    double t;
    if (!(hi > lo)) return 0;
    t = (v-lo)/(hi-lo);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    return (uint32_t)(t*levels) & levels;
}

int geomlib_spatial_sort2(const geomlib_vec2 *points, size_t count,
                          uint32_t *order)
{
    keyed *keys;
    double minx=DBL_MAX,miny=DBL_MAX,maxx=-DBL_MAX,maxy=-DBL_MAX;
    uint32_t side = 1, levels;
    size_t i;
    if (!points || !order || !count) return GEOMLIB_MESH_ERROR_INVALID;
    for (i = 0; i < count; ++i) {
        if (points[i].x<minx) minx=points[i].x;
        if (points[i].y<miny) miny=points[i].y;
        if (points[i].x>maxx) maxx=points[i].x;
        if (points[i].y>maxy) maxy=points[i].y;
    }
    while (side < 32768 && (double)side < maxx-minx+1.0) side *= 2;
    while (side < 32768 && (double)side < maxy-miny+1.0) side *= 2;
    levels = side-1;
    keys = (keyed *)malloc(count*sizeof(keyed));
    if (!keys) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i = 0; i < count; ++i) {
        uint32_t qx = quantize(points[i].x, minx, maxx, levels);
        uint32_t qy = quantize(points[i].y, miny, maxy, levels);
        keys[i].code = geomlib_hilbert2_encode(side, qx, qy);
        keys[i].id = (uint32_t)i;
    }
    qsort(keys, count, sizeof(keyed), cmp_keyed);
    for (i = 0; i < count; ++i) order[i] = keys[i].id;
    free(keys);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_spatial_sort3(const geomlib_vec3 *points, size_t count,
                          uint32_t *order)
{
    keyed *keys;
    double mnx=DBL_MAX,mny=DBL_MAX,mnz=DBL_MAX;
    double mxx=-DBL_MAX,mxy=-DBL_MAX,mxz=-DBL_MAX;
    const uint32_t levels = 65535;
    size_t i;
    if (!points || !order || !count) return GEOMLIB_MESH_ERROR_INVALID;
    for (i = 0; i < count; ++i) {
        if (points[i].x<mnx) mnx=points[i].x;
        if (points[i].y<mny) mny=points[i].y;
        if (points[i].z<mnz) mnz=points[i].z;
        if (points[i].x>mxx) mxx=points[i].x;
        if (points[i].y>mxy) mxy=points[i].y;
        if (points[i].z>mxz) mxz=points[i].z;
    }
    keys = (keyed *)malloc(count*sizeof(keyed));
    if (!keys) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i = 0; i < count; ++i) {
        uint32_t qx = quantize(points[i].x, mnx, mxx, levels);
        uint32_t qy = quantize(points[i].y, mny, mxy, levels);
        uint32_t qz = quantize(points[i].z, mnz, mxz, levels);
        keys[i].code = geomlib_morton3_encode(qx, qy, qz);
        keys[i].id = (uint32_t)i;
    }
    qsort(keys, count, sizeof(keyed), cmp_keyed);
    for (i = 0; i < count; ++i) order[i] = keys[i].id;
    free(keys);
    return GEOMLIB_MESH_SUCCESS;
}
