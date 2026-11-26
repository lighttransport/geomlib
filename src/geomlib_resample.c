/*
 * geomlib - Polyline resampling (see geomlib_resample.h)
 *
 * Clean-room cumulative-length interpolation. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_resample.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_resample_polyline(const geomlib_vec3 *points, size_t n,
                              size_t samples, geomlib_vec3 *out)
{
    double *cum;
    double total;
    size_t i, k;

    if (!points || !out || n<2 || samples<2) return GEOMLIB_MESH_ERROR_INVALID;
    cum=(double*)malloc(n*sizeof(double));
    if (!cum) return GEOMLIB_MESH_ERROR_MEMORY;
    cum[0]=0.0;
    for (i=1;i<n;++i)
        cum[i]=cum[i-1]+geomlib_vec3_length(geomlib_vec3_sub(points[i],points[i-1]));
    total=cum[n-1];
    if (!(total>1e-12)) {
        for (k=0;k<samples;++k) out[k]=points[0];
        free(cum);
        return GEOMLIB_MESH_SUCCESS;
    }
    k=0;
    for (i=0;i<samples;++i) {
        double target=total*i/(double)(samples-1);
        while (k+1<n-1 && cum[k+1]<target) ++k;
        {
            double seg=cum[k+1]-cum[k];
            double t=seg>1e-15 ? (target-cum[k])/seg : 0.0;
            if (t<0.0) t=0.0;
            if (t>1.0) t=1.0;
            out[i]=geomlib_vec3_add(points[k],
                geomlib_vec3_mul(geomlib_vec3_sub(points[k+1],points[k]),(float)t));
        }
    }
    out[0]=points[0];
    out[samples-1]=points[n-1];
    free(cum);
    return GEOMLIB_MESH_SUCCESS;
}
