/*
 * geomlib - Mean-value cage deformation (see geomlib_meancoords.h)
 *
 * Clean-room implementation of Floater's mean-value coordinates for planar
 * polygons. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_meancoords.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double cross2(geomlib_vec2 a, geomlib_vec2 b)
{
    return (double)a.x*b.y-(double)a.y*b.x;
}

static double dot2(geomlib_vec2 a, geomlib_vec2 b)
{
    return (double)a.x*b.x+(double)a.y*b.y;
}

/* Signed angle at q between rays to a and b, in (-pi, pi). */
static double signed_angle(geomlib_vec2 a, geomlib_vec2 q, geomlib_vec2 b)
{
    geomlib_vec2 u=geomlib_vec2_create(a.x-q.x,a.y-q.y);
    geomlib_vec2 v=geomlib_vec2_create(b.x-q.x,b.y-q.y);
    double cr=cross2(u,v), dt=dot2(u,v);
    double nu=hypot((double)u.x,u.y), nv=hypot((double)v.x,v.y);
    double ang;
    if (!(nu>1e-14 && nv>1e-14)) return 0.0;
    dt /= (nu*nv);
    if (dt>1.0) dt=1.0;
    if (dt<-1.0) dt=-1.0;
    ang=acos(dt);
    return cr>=0.0 ? ang : -ang;
}

int geomlib_mean_value_weights2d(const geomlib_vec2 *cage, size_t n,
                                 geomlib_vec2 q, double *weights)
{
    size_t i, j;
    double sum=0.0;
    if (!cage || !weights || n<3) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<n;++i) {
        if (fabs((double)cage[i].x-q.x) < 1e-9 &&
            fabs((double)cage[i].y-q.y) < 1e-9) {
            for (j=0;j<n;++j) weights[j]= j==i ? 1.0 : 0.0;
            return GEOMLIB_MESH_SUCCESS;
        }
    }
    for (i=0;i<n;++i) {
        geomlib_vec2 a,b,prev;
        double ra, wprev, wnext;
        j=(i+1)%n;
        prev=cage[(i+n-1)%n]; a=cage[i]; b=cage[j];
        ra=hypot((double)a.x-q.x,(double)a.y-q.y);
        if (!(ra>1e-14)) return GEOMLIB_MESH_ERROR_INVALID;
        /* alpha_{i-1} spans rays to prev and i; alpha_i spans i and next. */
        wprev=tan(0.5*fabs(signed_angle(prev,q,a)));
        wnext=tan(0.5*fabs(signed_angle(a,q,b)));
        weights[i]=(wprev+wnext)/ra;
        sum+=weights[i];
    }
    if (!(sum>1e-18)) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<n;++i) weights[i]/=sum;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_deform_mean_value2d(const geomlib_vec2 *cage,
                                const geomlib_vec2 *deformed, size_t n,
                                const geomlib_vec2 *points, size_t count,
                                geomlib_vec2 *out)
{
    double *w;
    size_t i, j;
    int rc;
    if (!cage || !deformed || !points || !out || n<3)
        return GEOMLIB_MESH_ERROR_INVALID;
    w=(double*)malloc(n*sizeof(double));
    if (!w) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i=0;i<count;++i) {
        double x=0.0,y=0.0;
        rc=geomlib_mean_value_weights2d(cage,n,points[i],w);
        if (rc!=GEOMLIB_MESH_SUCCESS) { free(w); return rc; }
        for (j=0;j<n;++j) {
            x+=w[j]*deformed[j].x;
            y+=w[j]*deformed[j].y;
        }
        out[i]=geomlib_vec2_create((float)x,(float)y);
    }
    free(w);
    return GEOMLIB_MESH_SUCCESS;
}
