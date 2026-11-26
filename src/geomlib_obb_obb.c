/*
 * geomlib - OBB / OBB overlap (see geomlib_obb_obb.h)
 *
 * Clean-room 15-axis separating-axis test. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_obb_obb.h"

#include <math.h>

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double o[3])
{
    o[0] = a[1] * b[2] - a[2] * b[1];
    o[1] = a[2] * b[0] - a[0] * b[2];
    o[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3(const double a[3])
{
    return sqrt(dot3(a, a));
}

static int load(const geomlib_obb *o, double w[3][3], double h[3], double c[3])
{
    int i;
    if (!o) return 0;
    for (i = 0; i < 3; ++i) {
        geomlib_vec3 ax = o->axes[i];
        double len = sqrt((double)ax.x * ax.x + (double)ax.y * ax.y +
                          (double)ax.z * ax.z);
        if (!(len > 1e-12)) return 0;
        w[i][0] = (double)ax.x / len;
        w[i][1] = (double)ax.y / len;
        w[i][2] = (double)ax.z / len;
        h[i] = fabs(i == 0 ? (double)o->half_extents.x
                    : (i == 1 ? (double)o->half_extents.y
                              : (double)o->half_extents.z));
        c[i] = i == 0 ? (double)o->center.x
             : i == 1 ? (double)o->center.y
                      : (double)o->center.z;
    }
    return 1;
}

int geomlib_obb_obb_overlap(const geomlib_obb *a, const geomlib_obb *b)
{
    double wa[3][3], wb[3][3], ha[3], hb[3], ca[3], cb[3], delta[3];
    double axis[3], ra, rb, d;
    int i, j;

    if (!load(a, wa, ha, ca) || !load(b, wb, hb, cb)) return 0;
    for (i = 0; i < 3; ++i) delta[i] = cb[i] - ca[i];

    for (j = 0; j < 15; ++j) {
        if (j < 3) {
            for (i = 0; i < 3; ++i) axis[i] = wa[j][i];
        } else if (j < 6) {
            for (i = 0; i < 3; ++i) axis[i] = wb[j - 3][i];
        } else {
            int p = (j - 6) / 3, q = (j - 6) % 3;
            cross3(wa[p], wb[q], axis);
        }
        if (norm3(axis) < 1e-12) continue;
        ra = 0.0;
        rb = 0.0;
        for (i = 0; i < 3; ++i) {
            ra += ha[i] * fabs(dot3(wa[i], axis));
            rb += hb[i] * fabs(dot3(wb[i], axis));
        }
        d = fabs(dot3(delta, axis));
        if (d > ra + rb + 1e-12) return 0;
    }
    return 1;
}
