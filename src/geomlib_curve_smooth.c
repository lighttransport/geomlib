/*
 * geomlib - Curve smoothing (see geomlib_curve_smooth.h)
 *
 * Clean-room Jacobi Laplacian flow with fixed endpoints. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_curve_smooth.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_smooth_polyline(geomlib_vec3 *points, size_t n,
                             int iterations, float lambda)
{
    geomlib_vec3 *next;
    size_t i;
    int it;

    if (!points || n<3 || iterations<0) return GEOMLIB_MESH_ERROR_INVALID;
    if (!(lambda>0.0f)) lambda=0.5f;
    if (lambda>1.0f) lambda=1.0f;
    next=(geomlib_vec3*)malloc(n*sizeof(geomlib_vec3));
    if (!next) return GEOMLIB_MESH_ERROR_MEMORY;
    for (it=0;it<iterations;++it) {
        next[0]=points[0];
        next[n-1]=points[n-1];
        for (i=1;i+1<n;++i) {
            geomlib_vec3 avg=geomlib_vec3_mul(
                geomlib_vec3_add(points[i-1],points[i+1]),0.5f);
            next[i]=geomlib_vec3_add(points[i],
                geomlib_vec3_mul(geomlib_vec3_sub(avg,points[i]),lambda));
        }
        memcpy(points,next,n*sizeof(geomlib_vec3));
    }
    free(next);
    return GEOMLIB_MESH_SUCCESS;
}
