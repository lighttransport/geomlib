/*
 * geomlib - Point in polygon (see geomlib_point_in_polygon.h)
 *
 * Clean-room ray crossing, winding number and boundary distance.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_point_in_polygon.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_point_in_polygon(const geomlib_vec2 *poly, size_t n, geomlib_vec2 q)
{
    size_t i;
    int inside=0;
    if (!poly||n<3) return 0;
    for (i=0;i<n;++i) {
        size_t j=(i+1)%n;
        if (((poly[i].y>q.y)!=(poly[j].y>q.y)) &&
            ((double)q.x < ((double)poly[j].x-poly[i].x)*((double)q.y-poly[i].y)/
                           ((double)poly[j].y-poly[i].y)+poly[i].x))
            inside=!inside;
    }
    return inside;
}

int geomlib_polygon_winding(const geomlib_vec2 *poly, size_t n, geomlib_vec2 q)
{
    size_t i;
    int wn=0;
    if (!poly||n<3) return 0;
    for (i=0;i<n;++i) {
        geomlib_vec2 a=poly[i],b=poly[(i+1)%n];
        double cr=((double)b.x-a.x)*((double)q.y-a.y)
                - ((double)b.y-a.y)*((double)q.x-a.x);
        if (a.y<=q.y) {
            if (b.y>q.y && cr>0.0) ++wn;
        } else {
            if (b.y<=q.y && cr<0.0) --wn;
        }
    }
    return wn;
}

static double seg_dist2(geomlib_vec2 q, geomlib_vec2 a, geomlib_vec2 b)
{
    double vx=(double)b.x-a.x,vy=(double)b.y-a.y;
    double wx=(double)q.x-a.x,wy=(double)q.y-a.y;
    double l2=vx*vx+vy*vy,t;
    geomlib_vec2 c;
    if (!(l2>1e-18)) return wx*wx+wy*wy;
    t=(wx*vx+wy*vy)/l2;
    if (t<0.0) t=0.0;
    if (t>1.0) t=1.0;
    c=geomlib_vec2_create((float)(a.x+t*(b.x-a.x)),
                          (float)(a.y+t*(b.y-a.y)));
    {
        double dx=(double)q.x-c.x,dy=(double)q.y-c.y;
        return dx*dx+dy*dy;
    }
}

float geomlib_polygon_signed_distance(const geomlib_vec2 *poly, size_t n,
                                      geomlib_vec2 q)
{
    size_t i;
    double best=DBL_MAX;
    if (!poly||n<3) return INFINITY;
    for (i=0;i<n;++i) {
        double d=seg_dist2(q,poly[i],poly[(i+1)%n]);
        if (d<best) best=d;
    }
    return (float)(geomlib_point_in_polygon(poly,n,q) ? -sqrt(best) : sqrt(best));
}
