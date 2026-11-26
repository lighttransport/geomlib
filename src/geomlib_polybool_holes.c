/*
 * geomlib - Polygon Boolean with holes (see geomlib_polybool_holes.h)
 *
 * Clean-room exact implementation: triangulate both polygons-with-holes and
 * combine convex triangle pieces (intersection via clipping, difference via
 * successive half-plane cuts). This preserves analytic area exactly.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polybool_holes.h"
#include "geomlib/geomlib_polytri_holes.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_polygon_with_holes geomlib_polygon_with_holes_create(void)
{
    geomlib_polygon_with_holes p;
    memset(&p,0,sizeof(p));
    return p;
}

void geomlib_polygon_with_holes_free(geomlib_polygon_with_holes *p)
{
    if (!p) return;
    free(p->vertices);
    free(p->offsets);
    memset(p,0,sizeof(*p));
}

int geomlib_polygon_with_holes_set(const geomlib_vec2 *outer, size_t no,
                                   const geomlib_vec2 *const *holes,
                                   const size_t *hole_sizes, size_t nh,
                                   geomlib_polygon_with_holes *out)
{
    geomlib_polygon_with_holes r=geomlib_polygon_with_holes_create();
    size_t total=no,i;
    if (!out||!outer||no<3||(nh&&(!holes||!hole_sizes))) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<nh;++i) total+=hole_sizes[i];
    r.offsets=(size_t*)calloc(nh+2,sizeof(size_t));
    r.vertices=(geomlib_vec2*)malloc(total*sizeof(geomlib_vec2));
    if (!r.offsets||!r.vertices) { geomlib_polygon_with_holes_free(&r);
        return GEOMLIB_MESH_ERROR_MEMORY; }
    memcpy(r.vertices,outer,no*sizeof(geomlib_vec2));
    r.offsets[0]=0; r.offsets[1]=no;
    for (i=0;i<nh;++i) {
        r.offsets[i+2]=r.offsets[i+1]+hole_sizes[i];
        memcpy(r.vertices+r.offsets[i+1],holes[i],hole_sizes[i]*sizeof(geomlib_vec2));
    }
    r.hole_count=nh;
    *out=r;
    return GEOMLIB_MESH_SUCCESS;
}

typedef struct {
    geomlib_vec2 *verts;
    size_t count,cap;
    size_t *offsets;
    size_t nloops,ocap;
} builder;

static int add_loop(builder *b, const geomlib_vec2 *p, size_t n)
{
    if (n<3) return 0;
    if (b->nloops+1>=b->ocap) {
        size_t nc=b->ocap?b->ocap*2:16;
        size_t *no=(size_t*)realloc(b->offsets,nc*sizeof(size_t));
        if (!no) return -1;
        b->offsets=no; b->ocap=nc;
    }
    if (b->count+n>b->cap) {
        size_t nc=b->cap?b->cap*2:64;
        geomlib_vec2 *nv;
        while (nc<b->count+n) nc*=2;
        nv=(geomlib_vec2*)realloc(b->verts,nc*sizeof(geomlib_vec2));
        if (!nv) return -1;
        b->verts=nv; b->cap=nc;
    }
    memcpy(b->verts+b->count,p,n*sizeof(geomlib_vec2));
    b->offsets[b->nloops]=b->count;
    b->count+=n;
    b->nloops++;
    b->offsets[b->nloops]=b->count;   /* keep the end sentinel valid */
    return 0;
}

static void builder_free(builder *b)
{
    free(b->verts); free(b->offsets);
    memset(b,0,sizeof(*b));
}

static void finish(builder *b, geomlib_poly_bool_result *out)
{
    if (b->nloops) {
        b->offsets[b->nloops]=b->count;
        out->loops=b->verts; b->verts=NULL;
        out->offsets=b->offsets; b->offsets=NULL;
        out->loop_count=b->nloops;
    }
    builder_free(b);
}

/* Clip a convex polygon by the half-plane cross(a,b,p) >= 0 (keep = sign). */
static size_t halfplane_clip(const geomlib_vec2 *in, size_t nin,
                             geomlib_vec2 a, geomlib_vec2 b, int keep_left,
                             geomlib_vec2 *out)
{
    size_t nout=0,i;
    for (i=0;i<nin;++i) {
        geomlib_vec2 s=in[i],e=in[(i+1)%nin];
        double ds=((double)b.x-a.x)*((double)s.y-a.y)
                - ((double)b.y-a.y)*((double)s.x-a.x);
        double de=((double)b.x-a.x)*((double)e.y-a.y)
                - ((double)b.y-a.y)*((double)e.x-a.x);
        int ins=keep_left ? ds>=-1e-10 : ds<=1e-10;
        int ine=keep_left ? de>=-1e-10 : de<=1e-10;
        if (ins) out[nout++]=s;
        if (ins!=ine) {
            double t=ds/(ds-de);
            out[nout++]=geomlib_vec2_create(
                (float)(s.x+t*(e.x-s.x)),(float)(s.y+t*(e.y-s.y)));
        }
    }
    return nout;
}

/* Convex intersection of two triangles. */
static size_t convex_intersect(const geomlib_vec2 *p, size_t np,
                               const geomlib_vec2 *q, size_t nq,
                               geomlib_vec2 *out)
{
    geomlib_vec2 tmp[64],cur[64];
    size_t nc=np,i;
    if (np>64||nq>64) return 0;
    memcpy(cur,p,np*sizeof(geomlib_vec2));
    for (i=0;i<nq && nc;++i) {
        nc=halfplane_clip(cur,nc,q[i],q[(i+1)%nq],1,tmp);
        memcpy(cur,tmp,nc*sizeof(geomlib_vec2));
    }
    memcpy(out,cur,nc*sizeof(geomlib_vec2));
    return nc;
}

/* Convex difference P \ Q as a set of convex pieces. */
static int convex_subtract(const geomlib_vec2 *p, size_t np,
                           const geomlib_vec2 *q, size_t nq, builder *b)
{
    geomlib_vec2 cur[64],piece[64];
    size_t nc=np,i;
    if (np>64||nq>64) return 0;
    memcpy(cur,p,np*sizeof(geomlib_vec2));
    for (i=0;i<nq && nc;++i) {
        size_t n=halfplane_clip(cur,nc,q[i],q[(i+1)%nq],0,piece);
        if (add_loop(b,piece,n)<0) return -1;
        nc=halfplane_clip(cur,nc,q[i],q[(i+1)%nq],1,piece);
        memcpy(cur,piece,nc*sizeof(geomlib_vec2));
    }
    return 0;
}

static int triangulate(const geomlib_polygon_with_holes *p,
                       geomlib_polytri_result *out)
{
    const geomlib_vec2 *holes[64];
    size_t sizes[64],i;
    if (p->hole_count>64) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<p->hole_count;++i) {
        holes[i]=p->vertices+p->offsets[i+1];
        sizes[i]=p->offsets[i+2]-p->offsets[i+1];
    }
    return geomlib_triangulate_with_holes(p->vertices,p->offsets[1],
        holes,sizes,p->hole_count,out);
}

int geomlib_polygon_boolean_holes(const geomlib_polygon_with_holes *a,
                                  const geomlib_polygon_with_holes *b,
                                  geomlib_poly_bool_op op,
                                  geomlib_poly_bool_result *out)
{
    geomlib_polytri_result ta=geomlib_polytri_result_create();
    geomlib_polytri_result tb=geomlib_polytri_result_create();
    builder bld;
    size_t i,j;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out||!a||!b) return rc;
    memset(&bld,0,sizeof(bld));
    rc=triangulate(a,&ta);
    if (rc!=GEOMLIB_MESH_SUCCESS) goto done;
    rc=triangulate(b,&tb);
    if (rc!=GEOMLIB_MESH_SUCCESS) goto done;

    if (op==GEOMLIB_POLY_INTERSECTION) {
        for (i=0;i<ta.triangle_count;++i)
            for (j=0;j<tb.triangle_count;++j) {
                geomlib_vec2 p[3],q[3],r[64];
                size_t k,n;
                for (k=0;k<3;++k) {
                    p[k]=ta.vertices[ta.triangles[i*3+k]];
                    q[k]=tb.vertices[tb.triangles[j*3+k]];
                }
                n=convex_intersect(p,3,q,3,r);
                if (add_loop(&bld,r,n)<0) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto done; }
            }
    } else if (op==GEOMLIB_POLY_DIFFERENCE) {
        for (i=0;i<ta.triangle_count;++i) {
            geomlib_vec2 p[3];
            size_t k;
            for (k=0;k<3;++k) p[k]=ta.vertices[ta.triangles[i*3+k]];
            {
                builder piece;
                memset(&piece,0,sizeof(piece));
                if (add_loop(&piece,p,3)<0) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto done; }
                /* Sequentially subtract each triangle of B from this piece. */
                for (j=0;j<tb.triangle_count;++j) {
                    geomlib_vec2 q[3];
                    size_t m;
                    for (m=0;m<3;++m) q[m]=tb.vertices[tb.triangles[j*3+m]];
                    {
                        builder next;
                        memset(&next,0,sizeof(next));
                        for (size_t l=0;l<piece.nloops;++l) {
                            size_t s=piece.offsets[l],e=piece.offsets[l+1];
                            if (convex_subtract(piece.verts+s,e-s,q,3,&next)<0) {
                                builder_free(&piece); builder_free(&next);
                                rc=GEOMLIB_MESH_ERROR_MEMORY; goto done;
                            }
                        }
                        builder_free(&piece);
                        piece=next;
                    }
                }
                for (size_t l=0;l<piece.nloops;++l)
                    if (add_loop(&bld,piece.verts+piece.offsets[l],
                                 piece.offsets[l+1]-piece.offsets[l])<0) {
                        builder_free(&piece);
                        rc=GEOMLIB_MESH_ERROR_MEMORY; goto done;
                    }
                builder_free(&piece);
            }
        }
    } else { /* Union: all of A plus B minus A */
        for (i=0;i<ta.triangle_count;++i) {
            geomlib_vec2 p[3];
            for (size_t k=0;k<3;++k) p[k]=ta.vertices[ta.triangles[i*3+k]];
            if (add_loop(&bld,p,3)<0) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto done; }
        }
        for (j=0;j<tb.triangle_count;++j) {
            geomlib_vec2 q[3];
            size_t k;
            for (k=0;k<3;++k) q[k]=tb.vertices[tb.triangles[j*3+k]];
            {
                builder piece,next;
                memset(&piece,0,sizeof(piece));
                memset(&next,0,sizeof(next));
                if (add_loop(&piece,q,3)<0) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto done; }
                for (i=0;i<ta.triangle_count;++i) {
                    geomlib_vec2 p[3];
                    for (k=0;k<3;++k) p[k]=ta.vertices[ta.triangles[i*3+k]];
                    memset(&next,0,sizeof(next));
                    for (size_t l=0;l<piece.nloops;++l) {
                        size_t s=piece.offsets[l],e=piece.offsets[l+1];
                        if (convex_subtract(piece.verts+s,e-s,p,3,&next)<0) {
                            builder_free(&piece); builder_free(&next);
                            rc=GEOMLIB_MESH_ERROR_MEMORY; goto done;
                        }
                    }
                    builder_free(&piece); piece=next;
                }
                for (size_t l=0;l<piece.nloops;++l)
                    if (add_loop(&bld,piece.verts+piece.offsets[l],
                                 piece.offsets[l+1]-piece.offsets[l])<0) {
                        builder_free(&piece);
                        rc=GEOMLIB_MESH_ERROR_MEMORY; goto done;
                    }
                builder_free(&piece);
            }
        }
    }
    finish(&bld,out);
    rc=GEOMLIB_MESH_SUCCESS;

done:
    if (rc!=GEOMLIB_MESH_SUCCESS) builder_free(&bld);
    geomlib_polytri_result_free(&ta);
    geomlib_polytri_result_free(&tb);
    return rc;
}
