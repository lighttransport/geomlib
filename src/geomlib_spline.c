/*
 * geomlib - Catmull-Rom resampling (see geomlib_spline.h)
 *
 * Clean-room Barry-Goldman pyramid evaluation. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_spline.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

size_t geomlib_catmull_output_count(size_t n, int closed,
                                    size_t samples_per_segment)
{
    size_t segs = closed ? n : (n ? n-1 : 0);
    return segs*samples_per_segment + 1;
}

static geomlib_vec3 lerp3(geomlib_vec3 a, geomlib_vec3 b, double t)
{
    return geomlib_vec3_add(geomlib_vec3_mul(a,(float)(1.0-t)),
                            geomlib_vec3_mul(b,(float)t));
}

/* Barry-Goldman pyramid with knots t0..t3, evaluating at t. */
static geomlib_vec3 evaluate(geomlib_vec3 p0, geomlib_vec3 p1,
                             geomlib_vec3 p2, geomlib_vec3 p3,
                             double t0, double t1, double t2, double t3,
                             double t)
{
    geomlib_vec3 a1, a2, a3, b1, b2, c;
    double d10=t1-t0, d21=t2-t1, d32=t3-t2;
    if (!(d10>1e-12)) { d10=1e-12; t1=t0+d10; }
    if (!(d21>1e-12)) { d21=1e-12; t2=t1+d21; }
    if (!(d32>1e-12)) { d32=1e-12; t3=t2+d32; }
    a1=lerp3(p0,p1,(t-t0)/d10);
    a2=lerp3(p1,p2,(t-t1)/d21);
    a3=lerp3(p2,p3,(t-t2)/d32);
    b1=lerp3(a1,a2,(t-t0)/(t2-t0));
    b2=lerp3(a2,a3,(t-t1)/(t3-t1));
    c=lerp3(b1,b2,(t-t1)/d21);
    return c;
}

int geomlib_resample_catmull_rom(const geomlib_vec3 *points, size_t n,
                                 int closed, size_t samples_per_segment,
                                 geomlib_catmull_type type,
                                 geomlib_vec3 *out)
{
    size_t segs, s, i, k=0;
    if (!points || !out || n<2 || !samples_per_segment)
        return GEOMLIB_MESH_ERROR_INVALID;
    segs = closed ? n : n-1;
    for (s=0;s<segs;++s) {
        uint32_t i0,i1,i2,i3;
        geomlib_vec3 p0,p1,p2,p3;
        double t0=0.0,t1=1.0,t2=2.0,t3=3.0;
        if (closed) {
            i0=(uint32_t)((s+n-1)%n);
            i1=(uint32_t)(s%n);
            i2=(uint32_t)((s+1)%n);
            i3=(uint32_t)((s+2)%n);
        } else {
            i0 = s? (uint32_t)(s-1) : 0;
            i1 = (uint32_t)s;
            i2 = (uint32_t)(s+1);
            i3 = (s+2<n) ? (uint32_t)(s+2) : (uint32_t)(n-1);
        }
        p0=points[i0]; p1=points[i1]; p2=points[i2]; p3=points[i3];
        if (type==GEOMLIB_CATMULL_CENTRIPETAL) {
            t1=t0+sqrt(geomlib_vec3_length(geomlib_vec3_sub(p1,p0)));
            t2=t1+sqrt(geomlib_vec3_length(geomlib_vec3_sub(p2,p1)));
            t3=t2+sqrt(geomlib_vec3_length(geomlib_vec3_sub(p3,p2)));
            /* Coincident controls must still give distinct knots. */
            if (!(t1>t0+1e-12)) t1=t0+1e-6;
            if (!(t2>t1+1e-12)) t2=t1+1e-6;
            if (!(t3>t2+1e-12)) t3=t2+1e-6;
        }
        for (i=0;i<samples_per_segment;++i) {
            double u=(double)i/samples_per_segment;
            double t=t1+(t2-t1)*u;
            out[k++]=evaluate(p0,p1,p2,p3,t0,t1,t2,t3,t);
        }
    }
    out[k]=closed ? points[0] : points[n-1];
    return GEOMLIB_MESH_SUCCESS;
}
