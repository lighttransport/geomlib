/*
 * geomlib - 2D harmonic cage coordinates (see geomlib_cage2d_harmonic.h)
 *
 * Clean-room angle-based weights with affine precision for convex cages.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_cage2d_harmonic.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double cross(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x-a.x)*((double)c.y-a.y)
         - ((double)b.y-a.y)*((double)c.x-a.x);
}

int geomlib_cage2d_harmonic_weights(const geomlib_vec2 *cage, size_t n,
                                    geomlib_vec2 q, double *weights)
{
    size_t i;
    double sum=0.0;
    if (!cage || !weights || n<3) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<n;++i) {
        geomlib_vec2 prev=cage[(i+n-1)%n],cur=cage[i],next=cage[(i+1)%n];
        geomlib_vec2 u=geomlib_vec2_create(prev.x-q.x,prev.y-q.y);
        geomlib_vec2 v=geomlib_vec2_create(cur.x-q.x,cur.y-q.y);
        geomlib_vec2 w=geomlib_vec2_create(next.x-q.x,next.y-q.y);
        double nu=sqrt(geomlib_vec2_dot(u,u)),nv=sqrt(geomlib_vec2_dot(v,v));
        double nw=sqrt(geomlib_vec2_dot(w,w));
        double ang1,ang2,c1,c2,w1,w2;
        if (!(nu>1e-12&&nv>1e-12&&nw>1e-12)) return GEOMLIB_MESH_ERROR_INVALID;
        c1=geomlib_vec2_dot(u,v)/(nu*nv);
        c2=geomlib_vec2_dot(v,w)/(nv*nw);
        if (c1>1.0) c1=1.0;
        if (c1<-1.0) c1=-1.0;
        if (c2>1.0) c2=1.0;
        if (c2<-1.0) c2=-1.0;
        ang1=acos(c1); ang2=acos(c2);
        /* Preserve sign from the cross product for non-convex configurations. */
        if (cross(q,prev,cur)<0.0) ang1=-ang1;
        if (cross(q,cur,next)<0.0) ang2=-ang2;
        w1=tan(0.5*ang1); w2=tan(0.5*ang2);
        weights[i]=(w1+w2)/nv;
        sum+=weights[i];
    }
    if (!(sum>1e-18)) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<n;++i) weights[i]/=sum;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_cage2d_harmonic_deform(const geomlib_vec2 *cage,
                                   const geomlib_vec2 *deformed, size_t n,
                                   const geomlib_vec2 *points, size_t count,
                                   geomlib_vec2 *out)
{
    double *w;
    size_t i,j;
    int rc;
    if (!cage||!deformed||!points||!out||n<3) return GEOMLIB_MESH_ERROR_INVALID;
    w=(double*)malloc(n*sizeof(double));
    if (!w) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i=0;i<count;++i) {
        double x=0.0,y=0.0;
        rc=geomlib_cage2d_harmonic_weights(cage,n,points[i],w);
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
