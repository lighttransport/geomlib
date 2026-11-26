/*
 * geomlib - 3D Poisson disk sampling (see geomlib_poisson3d.h)
 *
 * Clean-room Bridson dart throwing with deterministic PCG32 and a grid.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_poisson3d.h"
#include "geomlib/geomlib_random.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

size_t geomlib_poisson_disk_volume(const geomlib_poisson3d_params *params,
                                   geomlib_vec3 *out, size_t out_capacity)
{
    geomlib_poisson3d_params p;
    geomlib_pcg32 rng;
    double cell,r2;
    size_t nx,ny,nz,i,k,count=0;
    int32_t *grid=NULL;
    size_t *active=NULL;
    size_t na=0,ai=0;

    if (!params || !out || !out_capacity || !(params->radius>1e-6f))
        return 0;
    p=*params;
    if (p.candidates<1) p.candidates=30;
    if (!(p.box_max.x>p.box_min.x && p.box_max.y>p.box_min.y &&
          p.box_max.z>p.box_min.z)) return 0;
    if (!p.max_samples) p.max_samples=out_capacity;
    if (p.max_samples>out_capacity) p.max_samples=out_capacity;

    cell=p.radius/sqrt(3.0);
    nx=(size_t)ceil((p.box_max.x-p.box_min.x)/cell)+1;
    ny=(size_t)ceil((p.box_max.y-p.box_min.y)/cell)+1;
    nz=(size_t)ceil((p.box_max.z-p.box_min.z)/cell)+1;
    if (nx*ny*nz > (size_t)1<<26) return 0;
    grid=(int32_t*)malloc(nx*ny*nz*sizeof(int32_t));
    active=(size_t*)malloc(p.max_samples*sizeof(size_t));
    if (!grid||!active) { free(grid); free(active); return 0; }
    for (i=0;i<nx*ny*nz;++i) grid[i]=-1;
    r2=(double)p.radius*p.radius;

    p.seed = p.seed ? p.seed : 0x9e3779b9u;
    geomlib_pcg32_init(&rng,p.seed,(uint64_t)p.seed<<32);

    /* First sample at the box centre for determinism. */
    {
        geomlib_vec3 c=geomlib_vec3_create(
            (p.box_min.x+p.box_max.x)*0.5f,
            (p.box_min.y+p.box_max.y)*0.5f,
            (p.box_min.z+p.box_max.z)*0.5f);
        out[count]=c;
        grid[(size_t)((c.x-p.box_min.x)/cell)+
             nx*((size_t)((c.y-p.box_min.y)/cell)+
                 ny*(size_t)((c.z-p.box_min.z)/cell))]=(int32_t)count;
        active[na++]=count++;
    }

    while (ai<na && count<p.max_samples) {
        geomlib_vec3 c=out[active[ai]];
        int found=0;
        for (k=0;k<(size_t)p.candidates;++k) {
            double rad=p.radius*(1.0+geomlib_pcg32_nextf(&rng));
            double theta=2.0*M_PI*geomlib_pcg32_nextf(&rng);
            double z=2.0*geomlib_pcg32_nextf(&rng)-1.0;
            double s=sqrt(fmax(0.0,1.0-z*z));
            geomlib_vec3 q=geomlib_vec3_create(
                (float)(c.x+rad*s*cos(theta)),
                (float)(c.y+rad*s*sin(theta)),
                (float)(c.z+rad*z));
            size_t gx,gy,gz;
            int ok=1;
            int dx,dy,dz;
            if (q.x<p.box_min.x||q.x>p.box_max.x||
                q.y<p.box_min.y||q.y>p.box_max.y||
                q.z<p.box_min.z||q.z>p.box_max.z) continue;
            for (dx=-2;dx<=2 && ok;++dx)
                for (dy=-2;dy<=2 && ok;++dy)
                    for (dz=-2;dz<=2;++dz) {
                        long x=(long)((q.x-p.box_min.x)/cell)+dx;
                        long y=(long)((q.y-p.box_min.y)/cell)+dy;
                        long zz=(long)((q.z-p.box_min.z)/cell)+dz;
                        int32_t id;
                        if (x<0||y<0||zz<0||(size_t)x>=nx||(size_t)y>=ny||(size_t)zz>=nz)
                            continue;
                        id=grid[x+nx*(y+ny*zz)];
                        if (id>=0 && geomlib_vec3_length_squared(
                                geomlib_vec3_sub(out[id],q))<r2) { ok=0; break; }
                    }
            if (!ok) continue;
            gx=(size_t)((q.x-p.box_min.x)/cell);
            gy=(size_t)((q.y-p.box_min.y)/cell);
            gz=(size_t)((q.z-p.box_min.z)/cell);
            if (gx>=nx||gy>=ny||gz>=nz) continue;
            out[count]=q;
            grid[gx+nx*(gy+ny*gz)]=(int32_t)count;
            active[na++]=count++;
            found=1;
            break;
        }
        if (!found) { active[ai]=active[--na]; }
        else ai++;
    }
    free(grid); free(active);
    return count;
}
