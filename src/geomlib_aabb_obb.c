/*
 * geomlib - AABB / OBB overlap (see geomlib_aabb_obb.h)
 *
 * Clean-room 15-axis separating-axis test. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_aabb_obb.h"

#include <float.h>
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

int geomlib_aabb_obb_overlap(geomlib_vec3 amin, geomlib_vec3 amax,
                             const geomlib_obb *obb)
{
    double basis[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double w[3][3], hw[3], ha[3], delta[3], center[3];
    double axis[3], ra, rb, d;
    int i, j;

    if (!obb) return 0;
    for (i = 0; i < 3; ++i) {
        geomlib_vec3 ax = obb->axes[i];
        double len = sqrt((double)ax.x * ax.x + (double)ax.y * ax.y +
                          (double)ax.z * ax.z);
        if (!(len > 1e-12)) return 0;
        w[i][0] = (double)ax.x / len;
        w[i][1] = (double)ax.y / len;
        w[i][2] = (double)ax.z / len;
        hw[i] = fabs(i == 0 ? (double)obb->half_extents.x
                    : (i == 1 ? (double)obb->half_extents.y
                              : (double)obb->half_extents.z));
        ha[i] = 0.5 * fabs((double)(i == 0 ? amax.x - amin.x
                              : i == 1 ? amax.y - amin.y
                                       : amax.z - amin.z));
        center[i] = 0.5 * ((double)(i == 0 ? amin.x + amax.x
                              : i == 1 ? amin.y + amax.y
                                       : amin.z + amax.z));
        delta[i] = i == 0 ? (double)obb->center.x - center[i]
                 : i == 1 ? (double)obb->center.y - center[i]
                          : (double)obb->center.z - center[i];
    }

    for (j = 0; j < 15; ++j) {
        if (j < 3) {
            for (i = 0; i < 3; ++i) axis[i] = basis[j][i];
        } else if (j < 6) {
            for (i = 0; i < 3; ++i) axis[i] = w[j - 3][i];
        } else {
            int a = (j - 6) / 3, b = (j - 6) % 3;
            cross3(basis[a], w[b], axis);
        }
        if (norm3(axis) < 1e-12) continue;
        ra = 0.0;
        rb = 0.0;
        for (i = 0; i < 3; ++i) {
            ra += ha[i] * fabs((double)basis[i][0] * axis[0] +
                               (double)basis[i][1] * axis[1] +
                               (double)basis[i][2] * axis[2]);
            rb += hw[i] * fabs(dot3(w[i], axis));
        }
        d = fabs(dot3(delta, axis));
        if (d > ra + rb + 1e-12) return 0;
    }
    return 1;
}
