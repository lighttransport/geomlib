/*
 * geomlib - Geodesic path extraction (see geomlib_geodesic_path.h)
 *
 * Clean-room Dijkstra with predecessor reconstruction on the triangle-edge
 * graph. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_path.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t v; float d; } heap_item;

static void heap_push(heap_item *h, size_t *n, uint32_t v, float d)
{
    size_t i=(*n)++,p;
    h[i].v=v; h[i].d=d;
    while (i>0) {
        heap_item t;
        p=(i-1)/2;
        if (h[p].d<=h[i].d) break;
        t=h[p]; h[p]=h[i]; h[i]=t;
        i=p;
    }
}

static heap_item heap_pop(heap_item *h, size_t *n)
{
    heap_item top=h[0];
    --*n;
    if (*n) {
        size_t i=0;
        h[0]=h[*n];
        for (;;) {
            size_t l=2*i+1,r=l+1,m=i;
            heap_item t;
            if (l<*n && h[l].d<h[m].d) m=l;
            if (r<*n && h[r].d<h[m].d) m=r;
            if (m==i) break;
            t=h[m]; h[m]=h[i]; h[i]=t;
            i=m;
        }
    }
    return top;
}

geomlib_geodesic_path geomlib_geodesic_path_create(void)
{
    geomlib_geodesic_path p;
    memset(&p,0,sizeof(p));
    return p;
}

void geomlib_geodesic_path_free(geomlib_geodesic_path *p)
{
    if (!p) return;
    free(p->vertices);
    memset(p,0,sizeof(*p));
}

int geomlib_geodesic_path_between(const geomlib_trimesh *mesh,
                                  uint32_t source, uint32_t target,
                                  geomlib_geodesic_path *out)
{
    geomlib_geodesic_path r=geomlib_geodesic_path_create();
    uint32_t *adj_a=NULL,*prev=NULL;
    float *adj_w=NULL,*dist=NULL;
    size_t *off=NULL,*cursor=NULL,*adj_count=NULL;
    heap_item *heap=NULL;
    uint8_t *done=NULL;
    size_t n,i,j,k,e=0,hn=0,len=0;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->triangle_count || source>=mesh->vertex_count ||
        target>=mesh->vertex_count) return rc;
    n=mesh->vertex_count;
    if (source==target) {
        r.vertices=(uint32_t*)malloc(sizeof(uint32_t));
        if (!r.vertices) return GEOMLIB_MESH_ERROR_MEMORY;
        r.vertices[0]=source; r.count=1; r.length=0.0f;
        *out=r;
        return GEOMLIB_MESH_SUCCESS;
    }
    adj_count=(size_t*)calloc(n,sizeof(size_t));
    off=(size_t*)calloc(n+1,sizeof(size_t));
    dist=(float*)malloc(n*sizeof(float));
    prev=(uint32_t*)malloc(n*sizeof(uint32_t));
    done=(uint8_t*)calloc(n,sizeof(uint8_t));
    heap=(heap_item*)malloc((6*mesh->triangle_count+2)*sizeof(heap_item));
    if (!adj_count||!off||!dist||!prev||!done||!heap) goto cleanup;
    for (i=0;i<n;++i) { dist[i]=FLT_MAX; prev[i]=UINT32_MAX; }
    for (i=0;i<mesh->triangle_count;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        if (v[0]>=n||v[1]>=n||v[2]>=n) goto cleanup;
        for (k=0;k<3;++k) { adj_count[v[k]]++; adj_count[v[(k+1)%3]]++; }
    }
    for (i=0;i<n;++i) off[i+1]=off[i]+adj_count[i];
    e=off[n];
    adj_a=(uint32_t*)malloc(e*sizeof(uint32_t));
    adj_w=(float*)malloc(e*sizeof(float));
    cursor=(size_t*)malloc(n*sizeof(size_t));
    if (!adj_a||!adj_w||!cursor) goto cleanup;
    memcpy(cursor,off,n*sizeof(size_t));
    for (i=0;i<mesh->triangle_count;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        for (k=0;k<3;++k) {
            uint32_t a=v[k],b=v[(k+1)%3];
            float w=geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[a],
                                                         mesh->vertices[b]));
            size_t pa=cursor[a]++,pb=cursor[b]++;
            adj_a[pa]=b; adj_w[pa]=w;
            adj_a[pb]=a; adj_w[pb]=w;
        }
    }

    dist[source]=0.0f;
    heap_push(heap,&hn,source,0.0f);
    while (hn) {
        heap_item cur=heap_pop(heap,&hn);
        if (done[cur.v] || cur.d>dist[cur.v]) continue;
        done[cur.v]=1;
        if (cur.v==target) break;
        for (j=off[cur.v];j<off[cur.v+1];++j) {
            uint32_t nb=adj_a[j];
            float nd=cur.d+adj_w[j];
            if (nd<dist[nb]) {
                dist[nb]=nd; prev[nb]=cur.v;
                heap_push(heap,&hn,nb,nd);
            }
        }
    }
    if (!isfinite(dist[target])) goto cleanup;
    {
        uint32_t v=target;
        while (v!=UINT32_MAX) { len++; v=prev[v]; }
    }
    r.vertices=(uint32_t*)malloc(len*sizeof(uint32_t));
    if (!r.vertices) goto cleanup;
    {
        uint32_t v=target;
        for (j=len;j-- >0;) { r.vertices[j]=v; v=prev[v]; }
    }
    r.count=len;
    r.length=dist[target];
    *out=r;
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    free(adj_count); free(off); free(dist); free(prev); free(done); free(heap);
    free(adj_a); free(adj_w); free(cursor);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_geodesic_path_free(&r);
    return rc;
}
