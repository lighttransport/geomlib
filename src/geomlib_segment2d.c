/*
 * geomlib - 2D segment intersection (see geomlib_segment2d.h)
 *
 * Clean-room orientation-based classification. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_segment2d.h"

#include <float.h>
#include <math.h>
#include <string.h>

static double cross(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x-a.x)*((double)c.y-a.y)
         - ((double)b.y-a.y)*((double)c.x-a.x);
}

static int on_segment(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 p)
{
    return p.x>=fminf(a.x,b.x)-1e-9 && p.x<=fmaxf(a.x,b.x)+1e-9 &&
           p.y>=fminf(a.y,b.y)-1e-9 && p.y<=fmaxf(a.y,b.y)+1e-9;
}

static double dot2(geomlib_vec2 a, geomlib_vec2 b)
{
    return (double)a.x*b.x+(double)a.y*b.y;
}

geomlib_seg2d_result geomlib_segment2d_intersect(geomlib_vec2 a0, geomlib_vec2 a1,
                                                  geomlib_vec2 b0, geomlib_vec2 b1)
{
    geomlib_seg2d_result r;
    double d1=cross(a0,a1,b0),d2=cross(a0,a1,b1);
    double d3=cross(b0,b1,a0),d4=cross(b0,b1,a1);
    memset(&r,0,sizeof(r));
    r.kind=GEOMLIB_SEG2D_NONE;
    r.s=-1.0; r.t=-1.0;

    if (fabs(d1)<1e-12 && fabs(d2)<1e-12 &&
        fabs(d3)<1e-12 && fabs(d4)<1e-12) {
        /* Collinear: project endpoints onto A and keep the overlapping range. */
        geomlib_vec2 da=geomlib_vec2_sub(a1,a0);
        double len2=dot2(da,da);
        double sb=dot2(geomlib_vec2_sub(b0,a0),da)/len2;
        double se=dot2(geomlib_vec2_sub(b1,a0),da)/len2;
        double lo=fmax(0.0,fmin(sb,se)),hi=fmin(1.0,fmax(sb,se));
        if (lo>hi+1e-9) return r;
        r.kind=GEOMLIB_SEG2D_COLLINEAR;
        r.overlap[0]=geomlib_vec2_create((float)(a0.x+lo*da.x),
                                         (float)(a0.y+lo*da.y));
        r.overlap[1]=geomlib_vec2_create((float)(a0.x+hi*da.x),
                                         (float)(a0.y+hi*da.y));
        r.s=lo; r.t=sb+(se-sb)*lo;
        return r;
    }
    if (((d1>0&&d2<0)||(d1<0&&d2>0)) &&
        ((d3>0&&d4<0)||(d3<0&&d4>0))) {
        double s=d1/(d1-d2), t=d3/(d3-d4);
        r.kind=GEOMLIB_SEG2D_PROPER;
        r.point=geomlib_vec2_create((float)(a0.x+s*(a1.x-a0.x)),
                                    (float)(a0.y+s*(a1.y-a0.y)));
        r.s=s; r.t=t;
        return r;
    }
    if (fabs(d1)<1e-12 && on_segment(a0,a1,b0)) {
        r.kind=GEOMLIB_SEG2D_ENDPOINT; r.point=b0;
        r.s=dot2(geomlib_vec2_sub(b0,a0),geomlib_vec2_sub(a1,a0))/
            dot2(geomlib_vec2_sub(a1,a0),geomlib_vec2_sub(a1,a0));
        r.t=0.0;
    } else if (fabs(d2)<1e-12 && on_segment(a0,a1,b1)) {
        r.kind=GEOMLIB_SEG2D_ENDPOINT; r.point=b1;
        r.s=dot2(geomlib_vec2_sub(b1,a0),geomlib_vec2_sub(a1,a0))/
            dot2(geomlib_vec2_sub(a1,a0),geomlib_vec2_sub(a1,a0));
        r.t=1.0;
    } else if (fabs(d3)<1e-12 && on_segment(b0,b1,a0)) {
        r.kind=GEOMLIB_SEG2D_ENDPOINT; r.point=a0; r.s=0.0;
        r.t=dot2(geomlib_vec2_sub(a0,b0),geomlib_vec2_sub(b1,b0))/
            dot2(geomlib_vec2_sub(b1,b0),geomlib_vec2_sub(b1,b0));
    } else if (fabs(d4)<1e-12 && on_segment(b0,b1,a1)) {
        r.kind=GEOMLIB_SEG2D_ENDPOINT; r.point=a1; r.s=1.0;
        r.t=dot2(geomlib_vec2_sub(a1,b0),geomlib_vec2_sub(b1,b0))/
            dot2(geomlib_vec2_sub(b1,b0),geomlib_vec2_sub(b1,b0));
    }
    return r;
}
