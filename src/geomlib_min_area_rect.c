/*
 * geomlib - Minimum-area rectangle (see geomlib_min_area_rect.h)
 *
 * Clean-room rotating-calipers enumeration over convex hull edges.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_min_area_rect.h"
#include "geomlib/geomlib_hull2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_minimum_area_rectangle(const geomlib_vec2 *points, size_t n,
                                   geomlib_min_area_rect *out)
{
    geomlib_hull2d h=geomlib_hull2d_create();
    size_t i;
    double best=DBL_MAX;
    geomlib_min_area_rect result;

    if (!out || !points || n<2) return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_convex_hull2d(points,n,&h)!=GEOMLIB_MESH_SUCCESS)
        return GEOMLIB_MESH_ERROR_INVALID;
    memset(&result,0,sizeof(result));

    for (i=0;i<h.count;++i) {
        size_t j=(i+1)%h.count;
        geomlib_vec2 a=h.vertices[i],b=h.vertices[j];
        double ux=(double)b.x-a.x, uy=(double)b.y-a.y;
        double len=hypot(ux,uy);
        geomlib_vec2 u,v;
        double minu=DBL_MAX,maxu=-DBL_MAX,minv=DBL_MAX,maxv=-DBL_MAX;
        size_t k;
        if (!(len>1e-12)) continue;
        u=geomlib_vec2_create((float)(ux/len),(float)(uy/len));
        v=geomlib_vec2_create(-u.y,u.x);
        for (k=0;k<h.count;++k) {
            double dx=(double)h.vertices[k].x-a.x, dy=(double)h.vertices[k].y-a.y;
            double pu=dx*u.x+dy*u.y, pv=dx*v.x+dy*v.y;
            if (pu<minu) minu=pu;
            if (pu>maxu) maxu=pu;
            if (pv<minv) minv=pv;
            if (pv>maxv) maxv=pv;
        }
        {
            double wu=maxu-minu, wv=maxv-minv, ar=wu*wv;
            if (ar<best) {
                double cu=0.5*(minu+maxu), cv=0.5*(minv+maxv);
                best=ar;
                result.axes[0]=u; result.axes[1]=v;
                result.extents=geomlib_vec2_create((float)(0.5*wu),(float)(0.5*wv));
                result.center=geomlib_vec2_create(
                    (float)(a.x+cu*u.x+cv*v.x),
                    (float)(a.y+cu*u.y+cv*v.y));
                result.area=(float)ar;
            }
        }
    }
    geomlib_hull2d_free(&h);
    if (best==DBL_MAX) return GEOMLIB_MESH_ERROR_INVALID;
    *out=result;
    return GEOMLIB_MESH_SUCCESS;
}
