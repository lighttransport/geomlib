/*
 * geomlib - Edge-collapse coarsening (see geomlib_remesh_collapse.h)
 *
 * Clean-room implementation with link-condition and duplicate-face checks.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_remesh_collapse.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t a,b; int count; size_t t0,t1; } edge_rec;

static int cmp_edge(const void *pa, const void *pb)
{
    const edge_rec *x=(const edge_rec*)pa, *y=(const edge_rec*)pb;
    if (x->a<y->a) return -1;
    if (x->a>y->a) return 1;
    if (x->b<y->b) return -1;
    if (x->b>y->b) return 1;
    return 0;
}

static int cmp_tri(const void *pa, const void *pb)
{
    const uint32_t *x=(const uint32_t*)pa, *y=(const uint32_t*)pb;
    int i;
    for (i=0;i<3;++i) if (x[i]!=y[i]) return x[i]<y[i]?-1:1;
    return 0;
}

static int has_edge(const geomlib_trimesh *m, uint32_t a, uint32_t b)
{
    size_t i;
    for (i=0;i<m->triangle_count;++i) {
        uint32_t v[3]={m->indices[i*3],m->indices[i*3+1],m->indices[i*3+2]};
        int k;
        for (k=0;k<3;++k)
            if (((v[k]==a&&v[(k+1)%3]==b)||(v[k]==b&&v[(k+1)%3]==a))) return 1;
    }
    return 0;
}

static int is_neighbor(const geomlib_trimesh *m, uint32_t v, uint32_t x)
{
    size_t i;
    for (i=0;i<m->triangle_count;++i) {
        uint32_t t[3]={m->indices[i*3],m->indices[i*3+1],m->indices[i*3+2]};
        int k;
        for (k=0;k<3;++k)
            if ((t[k]==v && t[(k+1)%3]==x) || (t[k]==x && t[(k+1)%3]==v))
                return 1;
    }
    return 0;
}

static int link_ok(const geomlib_trimesh *m, uint32_t a, uint32_t b,
                   uint32_t c, uint32_t d)
{
    /* The link condition: common neighbours of a and b are exactly c and d. */
    size_t i;
    if (c==d || has_edge(m,c,d)) return 0;
    for (i=0;i<m->triangle_count;++i) {
        uint32_t v[3]={m->indices[i*3],m->indices[i*3+1],m->indices[i*3+2]};
        int k;
        for (k=0;k<3;++k) {
            uint32_t x=v[k], y=v[(k+1)%3];
            uint32_t other=0;
            if (x==a && y!=b && y!=a) other=y;
            else if (y==a && x!=b && x!=a) other=x;
            else continue;
            if (other==c || other==d) continue;
            if (is_neighbor(m,b,other)) return 0;
        }
    }
    return 1;
}

int geomlib_coarsen_edges(geomlib_trimesh *mesh, float min_edge_length,
                          int max_passes, int preserve_boundary)
{
    edge_rec *edges=NULL;
    int pass, collapses=0;
    size_t ne, i, j, unique;

    if (!mesh || !mesh->vertices || !mesh->indices || !mesh->triangle_count ||
        !(min_edge_length>0.0f)) return 0;
    if (max_passes<1) max_passes=1;
    edges=(edge_rec*)malloc(3*mesh->triangle_count*sizeof(edge_rec));
    if (!edges) return 0;

    for (pass=0;pass<max_passes;++pass) {
        uint32_t ca=UINT32_MAX,cb=UINT32_MAX;
        float best=FLT_MAX;
        edge_rec *ne_edges=(edge_rec*)realloc(edges,
            3*mesh->triangle_count*sizeof(edge_rec));
        if (!ne_edges) break;
        edges=ne_edges;
        ne=0;
        for (i=0;i<mesh->triangle_count;++i) {
            uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
            for (j=0;j<3;++j) {
                uint32_t a=v[j],b=v[(j+1)%3];
                float len=geomlib_vec3_length(geomlib_vec3_sub(
                    mesh->vertices[a],mesh->vertices[b]));
                edges[ne].a=a<b?a:b; edges[ne].b=a<b?b:a;
                edges[ne].count=1; edges[ne].t0=i; edges[ne++].t1=(size_t)-1;
                if (len>1e-6f && len<min_edge_length && len<best && a!=b) {
                    best=len; ca=edges[ne-1].a; cb=edges[ne-1].b;
                }
            }
        }
        if (ca==UINT32_MAX || best<=0.0f) break;
        qsort(edges,ne,sizeof(edge_rec),cmp_edge);
        unique=0;
        for (i=0;i<ne;) {
            j=i+1;
            while (j<ne && edges[j].a==edges[i].a && edges[j].b==edges[i].b) ++j;
            edges[unique]=edges[i];
            edges[unique].count=(int)(j-i);
            edges[unique].t1=(j-i==2)?edges[i+1].t0:(size_t)-1;
            unique++;
            i=j;
        }
        ne=unique;
        for (i=0;i<ne;++i)
            if (edges[i].a==ca && edges[i].b==cb) break;
        if (i==ne || edges[i].count!=2) break;
        if (preserve_boundary && edges[i].t1==(size_t)-1) break;
        {
            size_t ta=edges[i].t0,tb=edges[i].t1;
            uint32_t c=0,d=0,k;
            for (k=0;k<3;++k) {
                uint32_t x=mesh->indices[ta*3+k];
                if (x!=ca && x!=cb) c=x;
                x=mesh->indices[tb*3+k];
                if (x!=ca && x!=cb) d=x;
            }
            if (!link_ok(mesh,ca,cb,c,d)) break;
            mesh->vertices[ca]=geomlib_vec3_mul(
                geomlib_vec3_add(mesh->vertices[ca],mesh->vertices[cb]),0.5f);
            for (j=0;j<mesh->triangle_count*3;++j)
                if (mesh->indices[j]==cb) mesh->indices[j]=ca;
            /* Remove degenerate and duplicate triangles. */
            {
                size_t w=0;
                for (j=0;j<mesh->triangle_count;++j) {
                    uint32_t x=mesh->indices[j*3],y=mesh->indices[j*3+1],z=mesh->indices[j*3+2];
                    if (x==y||y==z||x==z) continue;
                    if (w!=j) memcpy(mesh->indices+w*3,mesh->indices+j*3,3*sizeof(uint32_t));
                    w++;
                }
                mesh->triangle_count=w;
            }
            qsort(mesh->indices,mesh->triangle_count,
                  sizeof(uint32_t)*3,cmp_tri);
            {
                size_t w=0;
                for (j=0;j<mesh->triangle_count;++j) {
                    if (w && memcmp(mesh->indices+w*3,mesh->indices+(w-1)*3,
                                    3*sizeof(uint32_t))==0) continue;
                    if (w!=j) memcpy(mesh->indices+w*3,mesh->indices+j*3,3*sizeof(uint32_t));
                    w++;
                }
                mesh->triangle_count=w;
            }
            collapses++;
        }
    }
    free(edges);

    /* Drop unreferenced vertices and remap indices. */
    {
        uint32_t *map=(uint32_t*)malloc(mesh->vertex_count*sizeof(uint32_t));
        geomlib_vec3 *nv;
        size_t i;
        if (map) {
            for (i=0;i<mesh->vertex_count;++i) map[i]=UINT32_MAX;
            for (i=0;i<mesh->triangle_count*3;++i)
                if (mesh->indices[i]<mesh->vertex_count) map[mesh->indices[i]]=0;
            {
                size_t w=0;
                nv=(geomlib_vec3*)malloc(mesh->vertex_count*sizeof(geomlib_vec3));
                if (nv) {
                    for (i=0;i<mesh->vertex_count;++i) {
                        if (map[i]==UINT32_MAX) continue;
                        map[i]=(uint32_t)w;
                        nv[w++]=mesh->vertices[i];
                    }
                    free(mesh->vertices);
                    mesh->vertices=nv;
                    mesh->vertex_count=w;
                    mesh->vertex_capacity=w;
                    for (i=0;i<mesh->triangle_count*3;++i)
                        mesh->indices[i]=map[mesh->indices[i]];
                }
            }
            free(map);
        }
    }
    return collapses;
}
