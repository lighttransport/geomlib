/*
 * geomlib - Dijkstra geodesic distances (see geomlib_geodesic_dijkstra.h)
 *
 * Clean-room binary-heap Dijkstra over the undirected triangle-edge graph.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_dijkstra.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t v; float w; } edge;
typedef struct { uint32_t v; float d; } heap_item;

static void heap_push(heap_item *h, size_t *n, uint32_t v, float d)
{
    size_t i = (*n)++;
    size_t p;
    h[i].v=v; h[i].d=d;
    while (i > 0) {
        heap_item t;
        p=(i-1)/2;
        if (h[p].d <= h[i].d) break;
        t=h[p]; h[p]=h[i]; h[i]=t;
        i=p;
    }
}

static heap_item heap_pop(heap_item *h, size_t *n)
{
    heap_item top = h[0];
    --*n;
    if (*n) {
        size_t i=0;
        h[0]=h[*n];
        for (;;) {
            size_t l=2*i+1, r=l+1, m=i;
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

int geomlib_geodesic_dijkstra(const geomlib_trimesh *mesh,
                              const uint32_t *sources, size_t source_count,
                              geomlib_geodesic_distances *out)
{
    edge *adj=NULL;
    size_t *offsets=NULL, *count=NULL, *cursor=NULL, ne=0, i, j, hn=0;
    heap_item *heap=NULL;
    uint8_t *done=NULL;
    size_t n;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->vertex_count || !sources || !source_count)
        return GEOMLIB_MESH_ERROR_INVALID;
    n=mesh->vertex_count;
    if (out->vertex_count != n || !out->distances)
        return GEOMLIB_MESH_ERROR_INVALID;

    count=(size_t*)calloc(n,sizeof(size_t));
    offsets=(size_t*)calloc(n+1,sizeof(size_t));
    if (!count||!offsets) goto fail;
    for (i=0;i<mesh->triangle_count;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        if (v[0]>=n||v[1]>=n||v[2]>=n) goto fail;
        for (j=0;j<3;++j) { count[v[j]]++; count[v[(j+1)%3]]++; }
    }
    for (i=0;i<n;++i) offsets[i+1]=offsets[i]+count[i];
    ne=offsets[n];
    adj=(edge*)malloc((ne?ne:1)*sizeof(edge));
    cursor=(size_t*)malloc(n*sizeof(size_t));
    heap=(heap_item*)malloc((6*mesh->triangle_count+2)*sizeof(heap_item));
    done=(uint8_t*)calloc(n,sizeof(uint8_t));
    if (!adj||!cursor||!heap||!done) goto fail;
    memcpy(cursor,offsets,n*sizeof(size_t));

    for (i=0;i<mesh->triangle_count;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        for (j=0;j<3;++j) {
            uint32_t a=v[j], b=v[(j+1)%3];
            float w=geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[a],
                                                         mesh->vertices[b]));
            adj[cursor[a]].v=b; adj[cursor[a]++].w=w;
            adj[cursor[b]].v=a; adj[cursor[b]++].w=w;
        }
    }

    for (i=0;i<n;++i) out->distances[i]=FLT_MAX;
    if (out->source_vertices) free(out->source_vertices);
    out->source_vertices=(uint32_t*)malloc(source_count*sizeof(uint32_t));
    if (!out->source_vertices) goto fail;
    for (i=0;i<source_count;++i) {
        out->source_vertices[i]=sources[i];
        if (sources[i]>=n) continue;
        out->distances[sources[i]]=0.0f;
        heap_push(heap,&hn,sources[i],0.0f);
    }
    out->source_count=source_count;

    while (hn) {
        heap_item cur=heap_pop(heap,&hn);
        if (done[cur.v] || cur.d > out->distances[cur.v]) continue;
        done[cur.v]=1;
        for (j=offsets[cur.v]; j<offsets[cur.v+1]; ++j) {
            uint32_t nb=adj[j].v;
            float nd=cur.d+adj[j].w;
            if (nd < out->distances[nb]) {
                out->distances[nb]=nd;
                heap_push(heap,&hn,nb,nd);
            }
        }
    }

    free(count); free(offsets); free(cursor); free(adj); free(heap); free(done);
    return GEOMLIB_MESH_SUCCESS;

fail:
    free(count); free(offsets); free(cursor); free(adj); free(heap); free(done);
    return GEOMLIB_MESH_ERROR_MEMORY;
}

int geomlib_geodesic_dijkstra_from_vertex(const geomlib_trimesh *mesh,
                                          uint32_t source,
                                          geomlib_geodesic_distances *out)
{
    return geomlib_geodesic_dijkstra(mesh,&source,1,out);
}
