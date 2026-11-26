/*
 * geomlib - Minimum spanning tree (see geomlib_mst.h)
 *
 * Clean-room Prim algorithm on the complete graph. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_mst.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_mst_result geomlib_mst_result_create(void)
{
    geomlib_mst_result r;
    memset(&r,0,sizeof(r));
    return r;
}

void geomlib_mst_result_free(geomlib_mst_result *r)
{
    if (!r) return;
    free(r->edges);
    memset(r,0,sizeof(*r));
}

static double dist2(const double *p,const double *q,int dim)
{
    double s=0.0;
    int k;
    for (k=0;k<dim;++k) {
        double d=p[k]-q[k];
        s+=d*d;
    }
    return s;
}

int geomlib_compute_mst(const double *points, size_t count, int dim,
                        geomlib_mst_result *out)
{
    geomlib_mst_result r=geomlib_mst_result_create();
    double *best=NULL;
    uint32_t *parent=NULL;
    uint8_t *used=NULL;
    size_t i,k;

    if (!out||!points||count<2||(dim!=2&&dim!=3))
        return GEOMLIB_MESH_ERROR_INVALID;
    best=(double*)malloc(count*sizeof(double));
    parent=(uint32_t*)malloc(count*sizeof(uint32_t));
    used=(uint8_t*)calloc(count,sizeof(uint8_t));
    r.edges=(uint32_t*)malloc((count-1)*2*sizeof(uint32_t));
    if (!best||!parent||!used||!r.edges) goto fail;
    for (i=0;i<count;++i) { best[i]=DBL_MAX; parent[i]=0; }
    best[0]=0.0;

    for (i=0;i<count;++i) {
        size_t v=(size_t)-1;
        double bd=DBL_MAX;
        for (k=0;k<count;++k)
            if (!used[k] && best[k]<bd) { bd=best[k]; v=k; }
        if (v==(size_t)-1) break;
        used[v]=1;
        if (v!=0) {
            r.edges[r.edge_count*2]=parent[v];
            r.edges[r.edge_count*2+1]=(uint32_t)v;
            r.edge_count++;
            r.total_weight+=sqrt(bd);
        }
        for (k=0;k<count;++k) {
            double d;
            if (used[k]) continue;
            d=dist2(points+v*dim,points+k*dim,dim);
            if (d<best[k]) { best[k]=d; parent[k]=(uint32_t)v; }
        }
    }
    if (r.edge_count!=count-1) goto fail;
    free(best); free(parent); free(used);
    *out=r;
    return GEOMLIB_MESH_SUCCESS;

fail:
    free(best); free(parent); free(used);
    geomlib_mst_result_free(&r);
    return GEOMLIB_MESH_ERROR_INVALID;
}
