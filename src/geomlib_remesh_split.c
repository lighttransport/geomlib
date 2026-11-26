/*
 * geomlib - Isotropic edge-split refinement (see geomlib_remesh_split.h)
 *
 * Clean-room implementation: split the longest offending edge into two
 * triangles (boundary) or four triangles (interior), rebuilding adjacency
 * after each split. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_remesh_split.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t a,b; size_t t0,t1; } edge_rec;

static int cmp_edge(const void *pa, const void *pb)
{
    const edge_rec *x=(const edge_rec*)pa, *y=(const edge_rec*)pb;
    if (x->a<y->a) return -1;
    if (x->a>y->a) return 1;
    if (x->b<y->b) return -1;
    if (x->b>y->b) return 1;
    return 0;
}

static geomlib_vec3 face_normal(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c)
{
    return geomlib_vec3_normalize(geomlib_vec3_cross(
        geomlib_vec3_sub(b,a),geomlib_vec3_sub(c,a)));
}

static void set_oriented(uint32_t *tri, uint32_t a, uint32_t b, uint32_t c,
                         geomlib_vec3 ref, const geomlib_vec3 *verts)
{
    geomlib_vec3 n=face_normal(verts[a],verts[b],verts[c]);
    if (geomlib_vec3_dot(n,ref)<0.0f) { uint32_t t=b; b=c; c=t; }
    tri[0]=a; tri[1]=b; tri[2]=c;
}

/* Append a vertex, returning its index or UINT32_MAX. */
static uint32_t append_vertex(geomlib_trimesh *m, geomlib_vec3 p)
{
    if (m->vertex_count == m->vertex_capacity) {
        size_t nc=m->vertex_capacity ? m->vertex_capacity*2 : 64;
        geomlib_vec3 *nv=(geomlib_vec3*)realloc(m->vertices,nc*sizeof(geomlib_vec3));
        if (!nv) return UINT32_MAX;
        m->vertices=nv; m->vertex_capacity=nc;
    }
    m->vertices[m->vertex_count]=p;
    return (uint32_t)m->vertex_count++;
}

/* Append a triangle, returning its index or SIZE_MAX. */
static size_t append_tri(geomlib_trimesh *m, uint32_t a, uint32_t b, uint32_t c)
{
    if (m->triangle_count == m->triangle_capacity) {
        size_t nc=m->triangle_capacity ? m->triangle_capacity*2 : 64;
        uint32_t *ni=(uint32_t*)realloc(m->indices,nc*3*sizeof(uint32_t));
        if (!ni) return (size_t)-1;
        m->indices=ni; m->triangle_capacity=nc;
    }
    m->indices[m->triangle_count*3]=a;
    m->indices[m->triangle_count*3+1]=b;
    m->indices[m->triangle_count*3+2]=c;
    return m->triangle_count++;
}

int geomlib_refine_edges(geomlib_trimesh *mesh, float max_edge_length,
                         int max_passes)
{
    edge_rec *edges=NULL;
    int pass, splits=0;
    size_t ne=0, i, j, unique;

    if (!mesh || !mesh->vertices || !mesh->indices || !mesh->triangle_count ||
        !(max_edge_length>0.0f)) return 0;
    if (max_passes<1) max_passes=1;
    edges=(edge_rec*)malloc(3*mesh->triangle_count*sizeof(edge_rec));
    if (!edges) return 0;

    for (pass=0;pass<max_passes;++pass) {
        uint32_t split_a=UINT32_MAX, split_b=UINT32_MAX;
        float split_len=0.0f;
        int did=0;
        {
            edge_rec *ne_edges=(edge_rec*)realloc(edges,
                3*mesh->triangle_count*sizeof(edge_rec));
            if (!ne_edges) break;
            edges=ne_edges;
        }
        ne=0;
        for (i=0;i<mesh->triangle_count;++i) {
            uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
            for (j=0;j<3;++j) {
                uint32_t a=v[j],b=v[(j+1)%3];
                float len;
                edges[ne].a = a<b?a:b;
                edges[ne].b = a<b?b:a;
                edges[ne].t0=i; edges[ne++].t1=(size_t)-1;
                len=geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[a],
                                                         mesh->vertices[b]));
                if (len>max_edge_length && len>split_len) {
                    split_len=len;
                    split_a=edges[ne-1].a; split_b=edges[ne-1].b;
                }
            }
        }
        if (split_a==UINT32_MAX) break;

        qsort(edges,ne,sizeof(edge_rec),cmp_edge);
        unique=0;
        for (i=0;i<ne;) {
            j=i+1;
            while (j<ne && edges[j].a==edges[i].a && edges[j].b==edges[i].b) ++j;
            edges[unique]=edges[i];
            edges[unique].t1 = (j-i==2) ? edges[i+1].t0 : (size_t)-1;
            unique++;
            i=j;
        }
        ne=unique;
        /* Locate the chosen edge after sorting by endpoints. */
        {
            uint32_t a=split_a,b=split_b;
            geomlib_vec3 mid;
            size_t t0,t1;
            uint32_t m,c=0,d=0;
            for (i=0;i<ne;++i)
                if (edges[i].a==a && edges[i].b==b) break;
            if (i==ne) break;
            t0=edges[i].t0; t1=edges[i].t1;
            for (j=0;j<3;++j) {
                uint32_t x=mesh->indices[t0*3+j];
                if (x!=a && x!=b) c=x;
                if (t1!=(size_t)-1) {
                    x=mesh->indices[t1*3+j];
                    if (x!=a && x!=b) d=x;
                }
            }
            mid=geomlib_vec3_mul(geomlib_vec3_add(mesh->vertices[a],mesh->vertices[b]),0.5f);
            m=append_vertex(mesh,mid);
            if (m==UINT32_MAX) break;
            {
                geomlib_vec3 ref=face_normal(mesh->vertices[mesh->indices[t0*3]],
                    mesh->vertices[mesh->indices[t0*3+1]],
                    mesh->vertices[mesh->indices[t0*3+2]]);
                set_oriented(mesh->indices+t0*3,a,m,c,ref,mesh->vertices);
                if (t1!=(size_t)-1) {
                    geomlib_vec3 ref2=face_normal(mesh->vertices[mesh->indices[t1*3]],
                        mesh->vertices[mesh->indices[t1*3+1]],
                        mesh->vertices[mesh->indices[t1*3+2]]);
                    set_oriented(mesh->indices+t1*3,b,m,d,ref2,mesh->vertices);
                    if (append_tri(mesh,m,b,c)==(size_t)-1) break;
                    if (append_tri(mesh,m,a,d)==(size_t)-1) break;
                    set_oriented(mesh->indices+(mesh->triangle_count-2)*3,
                                 m,b,c,ref,mesh->vertices);
                    set_oriented(mesh->indices+(mesh->triangle_count-1)*3,
                                 m,a,d,ref,mesh->vertices);
                } else {
                    if (append_tri(mesh,m,b,c)==(size_t)-1) break;
                    set_oriented(mesh->indices+(mesh->triangle_count-1)*3,
                                 m,b,c,ref,mesh->vertices);
                }
            }
            splits++; did=1;
        }
        if (!did) break;
    }
    free(edges);
    return splits;
}
