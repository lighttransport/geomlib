/*
 * geomlib - Edge-flip quality improvement (see geomlib_quality_flip.h)
 *
 * Clean-room implementation. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_quality_flip.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t a,b; size_t t0,t1; } edge_rec;

static int cmp_edge_rec(const void *pa, const void *pb)
{
    const edge_rec *x=(const edge_rec*)pa, *y=(const edge_rec*)pb;
    if (x->a<y->a) return -1;
    if (x->a>y->a) return 1;
    if (x->b<y->b) return -1;
    if (x->b>y->b) return 1;
    return 0;
}

static double corner_angle(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c)
{
    geomlib_vec3 u=geomlib_vec3_sub(b,a), v=geomlib_vec3_sub(c,a);
    float nu=geomlib_vec3_length(u), nv=geomlib_vec3_length(v);
    if (!(nu>1e-12f && nv>1e-12f)) return 0.0;
    {
        double cs=geomlib_vec3_dot(u,v)/(nu*nv);
        if (cs>1.0) cs=1.0;
        if (cs<-1.0) cs=-1.0;
        return acos(cs);
    }
}

static double tri_min_angle(const geomlib_trimesh *m, size_t t)
{
    uint32_t a=m->indices[t*3], b=m->indices[t*3+1], c=m->indices[t*3+2];
    geomlib_vec3 p=m->vertices[a], q=m->vertices[b], r=m->vertices[c];
    double x=corner_angle(p,q,r), y=corner_angle(q,r,p), z=corner_angle(r,p,q);
    return x<y ? (x<z?x:z) : (y<z?y:z);
}

float geomlib_mesh_min_angle(const geomlib_trimesh *mesh)
{
    size_t t;
    double mn=DBL_MAX;
    if (!mesh || !mesh->indices || !mesh->triangle_count) return 0.0f;
    for (t=0;t<mesh->triangle_count;++t) {
        double a=tri_min_angle(mesh,t);
        if (a<mn) mn=a;
    }
    return (float)mn;
}

static geomlib_vec3 face_normal(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c)
{
    return geomlib_vec3_normalize(geomlib_vec3_cross(
        geomlib_vec3_sub(b,a),geomlib_vec3_sub(c,a)));
}

/* Would flipping edge (a,b) with opposite vertices (c,d) keep orientation? */
static int convex_quad(const geomlib_trimesh *m, uint32_t a, uint32_t b,
                       uint32_t c, uint32_t d)
{
    /* Project onto the plane of the first triangle and require the two
     * diagonals (a,b) and (c,d) to cross in their interiors. */
    geomlib_vec3 pa=m->vertices[a], pb=m->vertices[b],
                 pc=m->vertices[c], pd=m->vertices[d];
    geomlib_vec3 n=face_normal(pa,pb,pc);
    geomlib_vec3 u=geomlib_vec3_normalize(geomlib_vec3_sub(pb,pa));
    geomlib_vec3 v=geomlib_vec3_cross(n,u);
    double ax=0.0,ay=0.0,bx,by,cx,cy,dx,dy;
    double o1,o2,o3,o4;
#define PX(p) (geomlib_vec3_dot(geomlib_vec3_sub((p),pa),u))
#define PY(p) (geomlib_vec3_dot(geomlib_vec3_sub((p),pa),v))
    bx=PX(pb); by=PY(pb);
    cx=PX(pc); cy=PY(pc);
    dx=PX(pd); dy=PY(pd);
    o1=ax*(by-cy)+bx*(cy-ay)+cx*(ay-by);
    o2=ax*(by-dy)+bx*(dy-ay)+dx*(ay-by);
    o3=cx*(dy-by)+dx*(by-cy)+bx*(cy-dy);
    o4=cx*(dy-ay)+dx*(ay-cy)+ax*(cy-dy);
#undef PX
#undef PY
    return o1*o2 < -1e-10 && o3*o4 < -1e-10;
}

int geomlib_improve_mesh_quality(geomlib_trimesh *mesh, int max_passes)
{
    edge_rec *edges=NULL;
    size_t ne=0, i, j, k;
    int pass, flips=0;

    if (!mesh || !mesh->vertices || !mesh->indices || !mesh->triangle_count)
        return 0;
    if (max_passes<1) max_passes=1;
    edges=(edge_rec*)malloc(3*mesh->triangle_count*sizeof(edge_rec));
    if (!edges) return 0;

    for (pass=0;pass<max_passes;++pass) {
        size_t unique=0;
        int changed=0;
        ne=0;
        for (i=0;i<mesh->triangle_count;++i) {
            uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
            for (k=0;k<3;++k) {
                uint32_t a=v[k], b=v[(k+1)%3];
                edges[ne].a = a<b?a:b;
                edges[ne].b = a<b?b:a;
                edges[ne].t0=i; edges[ne++].t1=(size_t)-1;
            }
        }
        qsort(edges,ne,sizeof(edge_rec),cmp_edge_rec);
        for (i=0;i<ne;) {
            j=i+1;
            while (j<ne && edges[j].a==edges[i].a && edges[j].b==edges[i].b) ++j;
            edges[unique]=edges[i];
            edges[unique].t1 = (j-i==2) ? edges[i+1].t0 : (size_t)-1;
            unique++;
            i=j;
        }
        ne=unique;

        for (i=0;i<ne;++i) {
            size_t ta=edges[i].t0, tb=edges[i].t1;
            uint32_t a=edges[i].a, b=edges[i].b, c=0, d=0;
            double before, after;
            if (ta>=mesh->triangle_count || tb>=mesh->triangle_count) continue;
            for (k=0;k<3;++k) {
                uint32_t x=mesh->indices[ta*3+k];
                if (x!=a && x!=b) c=x;
                x=mesh->indices[tb*3+k];
                if (x!=a && x!=b) d=x;
            }
            if (!c || !d || c==d) continue;
            /* Reject if edge (c,d) already exists. */
            for (j=0;j<ne;++j)
                if ((edges[j].a==c&&edges[j].b==d)||(edges[j].a==d&&edges[j].b==c))
                    break;
            if (j<ne) continue;
            if (!convex_quad(mesh,a,b,c,d)) continue;
            before=tri_min_angle(mesh,ta);
            { double t=tri_min_angle(mesh,tb); if (t<before) before=t; }
            {
                geomlib_vec3 p=mesh->vertices[a],q=mesh->vertices[d],r=mesh->vertices[c];
                geomlib_vec3 s=mesh->vertices[b];
                double x=corner_angle(p,q,r),y=corner_angle(q,r,p),z=corner_angle(r,p,q);
                after=x<y?(x<z?x:z):(y<z?y:z);
                x=corner_angle(s,r,q); y=corner_angle(r,q,s); z=corner_angle(q,s,r);
                { double t=x<y?(x<z?x:z):(y<z?y:z); if (t<after) after=t; }
            }
            if (after <= before+1e-4) continue;
            /* Flip the shared diagonal and orient the new faces consistently. */
            {
                uint32_t i0=mesh->indices[ta*3],i1=mesh->indices[ta*3+1],
                         i2=mesh->indices[ta*3+2];
                geomlib_vec3 ref=face_normal(mesh->vertices[i0],
                                             mesh->vertices[i1],mesh->vertices[i2]);
                geomlib_vec3 n=face_normal(mesh->vertices[a],mesh->vertices[d],
                                           mesh->vertices[c]);
                uint32_t x=d,y=c;
                if (geomlib_vec3_dot(ref,n)<0.0f) { x=c; y=d; }
                mesh->indices[ta*3]=a; mesh->indices[ta*3+1]=x; mesh->indices[ta*3+2]=y;
                n=face_normal(mesh->vertices[b],mesh->vertices[c],mesh->vertices[d]);
                x=c; y=d;
                if (geomlib_vec3_dot(ref,n)<0.0f) { x=d; y=c; }
                mesh->indices[tb*3]=b; mesh->indices[tb*3+1]=x; mesh->indices[tb*3+2]=y;
            }
            flips++; changed=1;
            /* The edge list is now stale; rebuild it on the next pass. */
            break;
        }
        if (!changed) break;
    }
    free(edges);
    return flips;
}
