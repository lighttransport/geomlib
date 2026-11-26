/*
 * geomlib - Region growing (see geomlib_region_growing.h)
 *
 * Clean-room BFS over face adjacency gated by dihedral angle.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_region_growing.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { uint32_t a,b; size_t f; } edge_rec;

static int cmp_edge(const void *pa, const void *pb)
{
    const edge_rec *x=(const edge_rec*)pa, *y=(const edge_rec*)pb;
    if (x->a<y->a) return -1;
    if (x->a>y->a) return 1;
    if (x->b<y->b) return -1;
    if (x->b>y->b) return 1;
    return 0;
}

static geomlib_vec3 face_normal(const geomlib_trimesh *m, size_t f)
{
    uint32_t a=m->indices[f*3],b=m->indices[f*3+1],c=m->indices[f*3+2];
    geomlib_vec3 n=geomlib_vec3_cross(geomlib_vec3_sub(m->vertices[b],m->vertices[a]),
                                      geomlib_vec3_sub(m->vertices[c],m->vertices[a]));
    float l=geomlib_vec3_length(n);
    return l>1e-12f ? geomlib_vec3_mul(n,1.0f/l) : geomlib_vec3_create(0,0,1);
}

int geomlib_mesh_region_growing(const geomlib_trimesh *mesh,
                                float dihedral_threshold_deg,
                                uint32_t *labels, size_t *region_count)
{
    edge_rec *edges=NULL;
    size_t *adj_head=NULL, *next=NULL;
    uint32_t *adj_f=NULL;
    geomlib_vec3 *normals=NULL;
    size_t ne=0, i, j, k, unique, nf;
    size_t regions=0;
    float cos_thresh;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !mesh->vertices || !mesh->indices || !mesh->triangle_count ||
        !labels) return rc;
    nf=mesh->triangle_count;
    if (!(dihedral_threshold_deg>0.0f)) dihedral_threshold_deg=30.0f;
    {
        double th=dihedral_threshold_deg*M_PI/180.0;
        cos_thresh=(float)cos(th);
    }
    for (i=0;i<nf;++i) labels[i]=UINT32_MAX;

    edges=(edge_rec*)malloc(3*nf*sizeof(edge_rec));
    normals=(geomlib_vec3*)malloc(nf*sizeof(geomlib_vec3));
    adj_head=(size_t*)malloc(nf*sizeof(size_t));
    next=(size_t*)malloc(6*nf*sizeof(size_t));
    adj_f=(uint32_t*)malloc(6*nf*sizeof(uint32_t));
    if (!edges||!normals||!adj_head||!next||!adj_f) goto cleanup;
    for (i=0;i<nf;++i) normals[i]=face_normal(mesh,i);
    for (i=0;i<nf;++i) adj_head[i]=(size_t)-1;

    for (i=0;i<nf;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        for (k=0;k<3;++k) {
            uint32_t a=v[k],b=v[(k+1)%3];
            edges[ne].a=a<b?a:b; edges[ne].b=a<b?b:a; edges[ne++].f=(uint32_t)i;
        }
    }
    qsort(edges,ne,sizeof(edge_rec),cmp_edge);
    unique=0;
    for (i=0;i<ne;) {
        j=i+1;
        while (j<ne && edges[j].a==edges[i].a && edges[j].b==edges[i].b) ++j;
        if (j-i==2) {
            size_t f0=edges[i].f, f1=edges[i+1].f;
            adj_f[unique]=(uint32_t)f1; next[unique]=adj_head[f0]; adj_head[f0]=unique++;
            adj_f[unique]=(uint32_t)f0; next[unique]=adj_head[f1]; adj_head[f1]=unique++;
        }
        i=j;
    }

    for (i=0;i<nf;++i) {
        uint32_t *queue;
        size_t head=0,tail=0;
        if (labels[i]!=UINT32_MAX) continue;
        labels[i]=(uint32_t)regions;
        queue=(uint32_t*)malloc(nf*sizeof(uint32_t));
        if (!queue) goto cleanup;
        queue[tail++]=(uint32_t)i;
        while (head<tail) {
            uint32_t f=queue[head++];
            size_t e;
            for (e=adj_head[f]; e!=(size_t)-1; e=next[e]) {
                uint32_t g=adj_f[e];
                float dot;
                if (labels[g]!=UINT32_MAX) continue;
                dot=geomlib_vec3_dot(normals[f],normals[g]);
                if (dot < cos_thresh) continue;
                labels[g]=(uint32_t)regions;
                queue[tail++]=(uint32_t)g;
            }
        }
        free(queue);
        regions++;
    }
    if (region_count) *region_count=regions;
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    free(edges); free(normals); free(adj_head); free(next); free(adj_f);
    return rc;
}
