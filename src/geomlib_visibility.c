/*
 * geomlib - Visibility graph shortest path (see geomlib_visibility.h)
 *
 * Clean-room visibility test + Dijkstra. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_visibility.h"
#include "geomlib/geomlib_segment2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_visibility_path geomlib_visibility_path_create(void)
{
    geomlib_visibility_path p;
    memset(&p,0,sizeof(p));
    return p;
}

void geomlib_visibility_path_free(geomlib_visibility_path *p)
{
    if (!p) return;
    free(p->path);
    memset(p,0,sizeof(*p));
}

static int visible(const geomlib_vec2 *pts, size_t np,
                   const geomlib_vec2 *const *obs, const size_t *sizes, size_t nh,
                   size_t i, size_t j)
{
    geomlib_vec2 a=pts[i],b=pts[j];
    size_t k,h;
    (void)np;
    for (h=0;h<nh;++h)
        for (k=0;k<sizes[h];++k) {
            geomlib_vec2 c=obs[h][k],d=obs[h][(k+1)%sizes[h]];
            geomlib_seg2d_result r=geomlib_segment2d_intersect(a,b,c,d);
            if (r.kind!=GEOMLIB_SEG2D_NONE) {
                /* Touching at a shared endpoint is allowed. */
                if (r.kind==GEOMLIB_SEG2D_ENDPOINT) {
                    if (fabs(r.point.x-a.x)<1e-8&&fabs(r.point.y-a.y)<1e-8) continue;
                    if (fabs(r.point.x-b.x)<1e-8&&fabs(r.point.y-b.y)<1e-8) continue;
                }
                /* Travelling along the obstacle edge itself is allowed. */
                if (r.kind==GEOMLIB_SEG2D_COLLINEAR) {
                    geomlib_seg2d_result ra=geomlib_segment2d_intersect(c,d,a,a);
                    geomlib_seg2d_result rb=geomlib_segment2d_intersect(c,d,b,b);
                    if (ra.kind!=GEOMLIB_SEG2D_NONE && rb.kind!=GEOMLIB_SEG2D_NONE)
                        continue;
                }
                return 0;
            }
        }
    return 1;
}

int geomlib_visibility_shortest_path(const geomlib_vec2 *start,
                                     const geomlib_vec2 *goal,
                                     const geomlib_vec2 *const *obstacles,
                                     const size_t *sizes, size_t nh,
                                     geomlib_visibility_path *out)
{
    geomlib_vec2 *pts=NULL;
    size_t np=2,off=2,i,j;
    double *dist=NULL,*w=NULL;
    int32_t *prev=NULL;
    uint8_t *done=NULL;
    geomlib_visibility_path r=geomlib_visibility_path_create();
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out||!start||!goal||(nh&&(!obstacles||!sizes)))
        return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<nh;++i) np+=sizes[i];
    pts=(geomlib_vec2*)malloc(np*sizeof(geomlib_vec2));
    if (!pts) return GEOMLIB_MESH_ERROR_MEMORY;
    pts[0]=*start; pts[1]=*goal;
    for (i=0;i<nh;++i) {
        memcpy(pts+off,obstacles[i],sizes[i]*sizeof(geomlib_vec2));
        off+=sizes[i];
    }
    dist=(double*)malloc(np*sizeof(double));
    w=(double*)malloc(np*np*sizeof(double));
    prev=(int32_t*)malloc(np*sizeof(int32_t));
    done=(uint8_t*)calloc(np,sizeof(uint8_t));
    if (!dist||!w||!prev||!done) goto cleanup;
    for (i=0;i<np;++i) {
        dist[i]=DBL_MAX; prev[i]=-1;
        for (j=0;j<np;++j) {
            w[i*np+j]=-1.0;
            if (i==j) { w[i*np+j]=0.0; continue; }
            if (visible(pts,np,obstacles,sizes,nh,i,j)) {
                double dx=(double)pts[i].x-pts[j].x;
                double dy=(double)pts[i].y-pts[j].y;
                w[i*np+j]=sqrt(dx*dx+dy*dy);
            }
        }
    }
    dist[0]=0.0;
    for (size_t it=0;it<np;++it) {
        size_t v=(size_t)-1;
        double bd=DBL_MAX;
        for (i=0;i<np;++i)
            if (!done[i]&&dist[i]<bd) { bd=dist[i]; v=i; }
        if (v==(size_t)-1) break;
        done[v]=1;
        for (i=0;i<np;++i)
            if (w[v*np+i]>0.0 && dist[v]+w[v*np+i]<dist[i]) {
                dist[i]=dist[v]+w[v*np+i];
                prev[i]=(int32_t)v;
            }
    }
    if (!(dist[1]<DBL_MAX)) goto cleanup;
    {
        size_t n=0;
        for (int32_t v=1;v!=-1;v=prev[v]) n++;
        r.path=(geomlib_vec2*)malloc(n*sizeof(geomlib_vec2));
        if (!r.path) goto cleanup;
        {
            size_t k=n;
            for (int32_t v=1;v!=-1;v=prev[v]) r.path[--k]=pts[v];
        }
        r.count=n;
        r.length=(float)dist[1];
    }
    *out=r;
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    free(pts); free(dist); free(w); free(prev); free(done);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_visibility_path_free(&r);
    return rc;
}
