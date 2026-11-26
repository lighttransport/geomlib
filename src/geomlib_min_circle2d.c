/*
 * geomlib - 2D minimum enclosing circle (see geomlib_min_circle2d.h)
 *
 * Clean-room Welzl incremental construction. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_min_circle2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { geomlib_vec2 c; float r; int valid; } circle;

static double dist(geomlib_vec2 a, geomlib_vec2 b)
{
    return hypot((double)a.x-b.x,(double)a.y-b.y);
}

static int contains(const circle *c, geomlib_vec2 p, double eps)
{
    return c->valid && dist(c->c,p) <= (double)c->r+eps;
}

static circle trivial(const geomlib_vec2 *p, size_t n)
{
    circle best,cur;
    unsigned mask;
    size_t i;
    best.valid=0; best.r=FLT_MAX; best.c=geomlib_vec2_create(0,0);
    for (mask=1;mask<(1u<<n);++mask) {
        geomlib_vec2 s[3];
        size_t m=0;
        int ok=1;
        for (i=0;i<n;++i) if (mask&(1u<<i)) s[m++]=p[i];
        cur.valid=0; cur.r=FLT_MAX;
        if (m==1) { cur.c=s[0]; cur.r=0.0f; cur.valid=1; }
        else if (m==2) {
            cur.c=geomlib_vec2_create((s[0].x+s[1].x)*0.5f,
                                      (s[0].y+s[1].y)*0.5f);
            cur.r=(float)(0.5*dist(s[0],s[1]));
            cur.valid=1;
        } else if (m==3) {
            double ax=s[0].x,ay=s[0].y,bx=s[1].x,by=s[1].y;
            double cx=s[2].x,cy=s[2].y;
            double d=2.0*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
            if (fabs(d)>1e-14) {
                double a2=ax*ax+ay*ay,b2=bx*bx+by*by,c2=cx*cx+cy*cy;
                cur.c=geomlib_vec2_create(
                    (float)((a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d),
                    (float)((a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d));
                cur.r=(float)dist(cur.c,s[0]);
                cur.valid=1;
            }
        }
        if (!cur.valid) continue;
        for (i=0;i<n;++i) if (!contains(&cur,p[i],1e-7)) { ok=0; break; }
        if (ok && cur.r<best.r) best=cur;
    }
    return best;
}

static circle welzl(const geomlib_vec2 *p, size_t n,
                    const geomlib_vec2 *r, size_t nr)
{
    circle d;
    size_t i;
    if (n==0 || nr==3) return trivial(r,nr);
    d=welzl(p,n-1,r,nr);
    if (contains(&d,p[n-1],1e-7)) return d;
    {
        geomlib_vec2 next[3];
        for (i=0;i<nr;++i) next[i]=r[i];
        next[nr]=p[n-1];
        return welzl(p,n-1,next,nr+1);
    }
}

int geomlib_minimum_enclosing_circle(const geomlib_vec2 *points, size_t count,
                                     geomlib_min_circle2d *out)
{
    circle c;
    if (!out || !points || !count) return GEOMLIB_MESH_ERROR_INVALID;
    if (count>512) return GEOMLIB_MESH_ERROR_INVALID;
    c=welzl(points,count,NULL,0);
    if (!c.valid) return GEOMLIB_MESH_ERROR_INVALID;
    out->center=c.c; out->radius=c.r;
    return GEOMLIB_MESH_SUCCESS;
}
