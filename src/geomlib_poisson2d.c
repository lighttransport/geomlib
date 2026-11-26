/*
 * geomlib - 2D Poisson disk in a polygon (see geomlib_poisson2d.h)
 *
 * Clean-room Bridson dart throwing with grid acceleration and polygon test.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_poisson2d.h"
#include "geomlib/geomlib_point_in_polygon.h"
#include "geomlib/geomlib_random.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

size_t geomlib_poisson_disk_polygon(const geomlib_vec2 *polygon, size_t np,
                                    const geomlib_poisson2d_params *params,
                                    geomlib_vec2 *out, size_t out_capacity)
{
    geomlib_poisson2d_params p;
    geomlib_pcg32 rng;
    double cell,r2,bminx,bminy,bmaxx,bmaxy;
    size_t nx,ny,i,k,count=0;
    int32_t *grid=NULL;
    size_t *active=NULL;
    size_t na=0,ai=0;

    if (!params||!polygon||np<3||!out||!out_capacity||!(params->radius>1e-6f))
        return 0;
    p=*params;
    if (p.candidates<1) p.candidates=30;
    if (!p.max_samples) p.max_samples=out_capacity;
    if (p.max_samples>out_capacity) p.max_samples=out_capacity;
    bminx=bminy=DBL_MAX; bmaxx=bmaxy=-DBL_MAX;
    for (i=0;i<np;++i) {
        if (polygon[i].x<bminx) bminx=polygon[i].x;
        if (polygon[i].y<bminy) bminy=polygon[i].y;
        if (polygon[i].x>bmaxx) bmaxx=polygon[i].x;
        if (polygon[i].y>bmaxy) bmaxy=polygon[i].y;
    }
    cell=p.radius/sqrt(2.0);
    nx=(size_t)ceil((bmaxx-bminx)/cell)+1;
    ny=(size_t)ceil((bmaxy-bminy)/cell)+1;
    if (nx*ny>(size_t)1<<26) return 0;
    grid=(int32_t*)malloc(nx*ny*sizeof(int32_t));
    active=(size_t*)malloc(p.max_samples*sizeof(size_t));
    if (!grid||!active) { free(grid); free(active); return 0; }
    for (i=0;i<nx*ny;++i) grid[i]=-1;
    r2=(double)p.radius*p.radius;
    p.seed=p.seed?p.seed:0x9e3779b9u;
    geomlib_pcg32_init(&rng,p.seed,(uint64_t)p.seed<<32);

    /* Seed at the first polygon vertex that lies inside. */
    {
        int seeded=0;
        for (i=0;i<np;++i) {
            if (geomlib_point_in_polygon(polygon,np,polygon[i])) {
                /* nudged inward toward the centroid */
                double cx=0.0,cy=0.0;
                size_t j;
                for (j=0;j<np;++j) { cx+=polygon[j].x; cy+=polygon[j].y; }
                cx/=np; cy/=np;
                out[0]=geomlib_vec2_create(
                    (float)(polygon[i].x+0.1*(cx-polygon[i].x)),
                    (float)(polygon[i].y+0.1*(cy-polygon[i].y)));
                if (!geomlib_point_in_polygon(polygon,np,out[0])) {
                    out[0]=polygon[i];
                }
                grid[(size_t)((out[0].x-bminx)/cell)+
                     nx*(size_t)((out[0].y-bminy)/cell)]=(int32_t)count;
                active[na++]=count++;
                seeded=1;
                break;
            }
        }
        if (!seeded) { free(grid); free(active); return 0; }
    }

    while (ai<na && count<p.max_samples) {
        geomlib_vec2 c=out[active[ai]];
        int found=0;
        for (k=0;k<(size_t)p.candidates;++k) {
            double rad=p.radius*(1.0+geomlib_pcg32_nextf(&rng));
            double ang=2.0*M_PI*geomlib_pcg32_nextf(&rng);
            geomlib_vec2 q=geomlib_vec2_create(
                (float)(c.x+rad*cos(ang)),(float)(c.y+rad*sin(ang)));
            size_t gx,gy;
            int ok=1,dx,dy;
            if (!geomlib_point_in_polygon(polygon,np,q)) continue;
            for (dx=-2;dx<=2&&ok;++dx)
                for (dy=-2;dy<=2;++dy) {
                    long x=(long)((q.x-bminx)/cell)+dx;
                    long y=(long)((q.y-bminy)/cell)+dy;
                    int32_t id;
                    if (x<0||y<0||(size_t)x>=nx||(size_t)y>=ny) continue;
                    id=grid[x+nx*y];
                    if (id>=0 && geomlib_vec2_length_squared(
                            geomlib_vec2_sub(out[id],q))<r2) { ok=0; break; }
                }
            if (!ok) continue;
            gx=(size_t)((q.x-bminx)/cell); gy=(size_t)((q.y-bminy)/cell);
            if (gx>=nx||gy>=ny) continue;
            out[count]=q;
            grid[gx+nx*gy]=(int32_t)count;
            active[na++]=count++;
            found=1;
            break;
        }
        if (!found) active[ai]=active[--na];
    }
    free(grid); free(active);
    return count;
}
