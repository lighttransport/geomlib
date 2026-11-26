/*
 * geomlib - Largest inscribed circle (see geomlib_inscribed_circle.h)
 *
 * Clean-room sampled-boundary Voronoi evaluation. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_inscribed_circle.h"
#include "geomlib/geomlib_voronoi2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double seg_dist(geomlib_vec2 q, geomlib_vec2 a, geomlib_vec2 b)
{
    double vx=(double)b.x-a.x, vy=(double)b.y-a.y;
    double wx=(double)q.x-a.x, wy=(double)q.y-a.y;
    double l2=vx*vx+vy*vy, t;
    geomlib_vec2 c;
    if (!(l2>1e-18)) return hypot((double)q.x-a.x,(double)q.y-a.y);
    t=(wx*vx+wy*vy)/l2;
    if (t<0.0) t=0.0;
    if (t>1.0) t=1.0;
    c=geomlib_vec2_create((float)(a.x+t*(b.x-a.x)),(float)(a.y+t*(b.y-a.y)));
    return hypot((double)q.x-c.x,(double)q.y-c.y);
}

static double boundary_distance(geomlib_vec2 q, const geomlib_vec2 *p, size_t n)
{
    size_t k;
    double best=DBL_MAX;
    for (k=0;k<n;++k) {
        double d=seg_dist(q,p[k],p[(k+1)%n]);
        if (d<best) best=d;
    }
    return best;
}

static int inside(geomlib_vec2 q, const geomlib_vec2 *p, size_t n)
{
    size_t k;
    int in=0;
    for (k=0;k<n;++k) {
        size_t j=(k+1)%n;
        if (((p[k].y>q.y)!=(p[j].y>q.y)) &&
            ((double)q.x < ((double)p[j].x-p[k].x)*((double)q.y-p[k].y)/
                           ((double)p[j].y-p[k].y)+p[k].x))
            in=!in;
    }
    return in;
}

int geomlib_largest_inscribed_circle(const geomlib_vec2 *polygon, size_t n,
                                     float spacing,
                                     geomlib_inscribed_circle *out)
{
    geomlib_vec2 *samples=NULL;
    geomlib_voronoi2d vor=geomlib_voronoi2d_create();
    size_t np=0,k;
    double minx=DBL_MAX,miny=DBL_MAX,maxx=-DBL_MAX,maxy=-DBL_MAX,diag;
    double best=-1.0;
    geomlib_vec2 bestc={0,0};
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !polygon || n<3) return rc;
    for (k=0;k<n;++k) {
        if (polygon[k].x<minx) minx=polygon[k].x;
        if (polygon[k].y<miny) miny=polygon[k].y;
        if (polygon[k].x>maxx) maxx=polygon[k].x;
        if (polygon[k].y>maxy) maxy=polygon[k].y;
    }
    diag=hypot(maxx-minx,maxy-miny);
    if (!(diag>1e-9)) return rc;
    if (!(spacing>0.0f)) spacing=(float)(0.02*diag);
    if (spacing < diag*1e-4f) spacing=(float)(diag*1e-4f);

    for (k=0;k<n;++k) {
        geomlib_vec2 a=polygon[k],b=polygon[(k+1)%n];
        double len=hypot((double)b.x-a.x,(double)b.y-a.y);
        int parts=(int)ceil(len/spacing),t;
        if (parts<1) parts=1;
        samples=(geomlib_vec2*)realloc(samples,(np+(size_t)parts+1)*sizeof(geomlib_vec2));
        if (!samples) goto cleanup;
        for (t=0;t<parts;++t) {
            double f=(double)t/parts;
            samples[np++]=geomlib_vec2_create((float)(a.x+f*(b.x-a.x)),
                                              (float)(a.y+f*(b.y-a.y)));
        }
    }
    {
        geomlib_vec2 lo=geomlib_vec2_create((float)minx,(float)miny);
        geomlib_vec2 hi=geomlib_vec2_create((float)maxx,(float)maxy);
        if (geomlib_compute_voronoi2d(samples,np,lo,hi,&vor)!=GEOMLIB_MESH_SUCCESS)
            goto cleanup;
    }
    for (k=0;k<vor.vertex_count;++k) {
        geomlib_vec2 q=vor.cells[k];
        double d;
        if (!inside(q,polygon,n)) continue;
        d=boundary_distance(q,polygon,n);
        if (d>best) { best=d; bestc=q; }
    }
    if (best<=0.0) goto cleanup;
    out->center=bestc;
    out->radius=(float)best;
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    free(samples);
    geomlib_voronoi2d_free(&vor);
    return rc;
}
