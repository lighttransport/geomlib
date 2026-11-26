/*
 * geomlib - Polygon rasterization (see geomlib_raster2d.h)
 *
 * Clean-room cell-centre sampling with point-in-polygon and signed distance.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_raster2d.h"
#include "geomlib/geomlib_point_in_polygon.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_raster2d geomlib_raster2d_create(void)
{
    geomlib_raster2d r;
    memset(&r,0,sizeof(r));
    return r;
}

void geomlib_raster2d_free(geomlib_raster2d *r)
{
    if (!r) return;
    free(r->mask);
    free(r->distance);
    memset(r,0,sizeof(*r));
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
    c=geomlib_vec2_create((float)(a.x+t*(b.x-a.x)),(float)(a.y+t*(b.y-a.y)));
    {
        double dx=(double)q.x-c.x,dy=(double)q.y-c.y;
        return dx*dx+dy*dy;
    }
}

int geomlib_rasterize_polygon2d(const geomlib_vec2 *polygon, size_t n,
                                uint32_t nx, uint32_t ny, float margin,
                                int with_distance, geomlib_raster2d *out)
{
    geomlib_raster2d r=geomlib_raster2d_create();
    double bminx=DBL_MAX,bminy=DBL_MAX,bmaxx=-DBL_MAX,bmaxy=-DBL_MAX;
    size_t i,x,y;

    if (!out||!polygon||n<3||!nx||!ny) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<n;++i) {
        if (polygon[i].x<bminx) bminx=polygon[i].x;
        if (polygon[i].y<bminy) bminy=polygon[i].y;
        if (polygon[i].x>bmaxx) bmaxx=polygon[i].x;
        if (polygon[i].y>bmaxy) bmaxy=polygon[i].y;
    }
    if (!(bmaxx>bminx)||!(bmaxy>bminy)) return GEOMLIB_MESH_ERROR_INVALID;
    if (!(margin>=0.0f)) margin=0.05f;
    r.nx=nx; r.ny=ny;
    r.cell_size=(float)fmax((bmaxx-bminx)*(1.0+2.0*margin)/nx,
                            (bmaxy-bminy)*(1.0+2.0*margin)/ny);
    r.origin=geomlib_vec2_create(
        (float)(0.5*(bminx+bmaxx)-0.5*nx*r.cell_size),
        (float)(0.5*(bminy+bmaxy)-0.5*ny*r.cell_size));
    r.mask=(uint8_t*)calloc((size_t)nx*ny,1);
    if (!r.mask) return GEOMLIB_MESH_ERROR_MEMORY;
    if (with_distance) {
        r.distance=(float*)malloc((size_t)nx*ny*sizeof(float));
        if (!r.distance) { geomlib_raster2d_free(&r); return GEOMLIB_MESH_ERROR_MEMORY; }
    }
    for (y=0;y<ny;++y)
        for (x=0;x<nx;++x) {
            geomlib_vec2 q=geomlib_vec2_create(
                r.origin.x+((float)x+0.5f)*r.cell_size,
                r.origin.y+((float)y+0.5f)*r.cell_size);
            size_t id=x+nx*y;
            int inside=geomlib_point_in_polygon(polygon,n,q);
            r.mask[id]=(uint8_t)inside;
            if (r.distance) {
                double best=DBL_MAX;
                for (i=0;i<n;++i) {
                    double d=seg_dist2(q,polygon[i],polygon[(i+1)%n]);
                    if (d<best) best=d;
                }
                r.distance[id]=(float)(inside ? -sqrt(best) : sqrt(best));
            }
        }
    *out=r;
    return GEOMLIB_MESH_SUCCESS;
}
