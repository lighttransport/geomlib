/*
 * geomlib - Broad-phase overlap queries (see geomlib_mesh_overlap.h)
 *
 * Clean-room BVH-style box tree with pair enumeration. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_mesh_overlap.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { geomlib_aabb box; uint32_t lo,hi; int left,right; } ov_node;

typedef struct {
    ov_node *nodes;
    size_t count,cap;
} ov_tree;

static int boxes_overlap(geomlib_aabb a, geomlib_aabb b)
{
    return !(a.max.x < b.min.x || b.max.x < a.min.x ||
             a.max.y < b.min.y || b.max.y < a.min.y ||
             a.max.z < b.min.z || b.max.z < a.min.z);
}

static geomlib_aabb tri_box(const geomlib_trimesh *m, size_t t)
{
    geomlib_vec3 p[3];
    size_t i;
    for (i=0;i<3;++i) {
        uint32_t id=m->indices[t*3+i];
        p[i]= id<m->vertex_count ? m->vertices[id] : geomlib_vec3_create(0,0,0);
    }
    {
        geomlib_vec3 lo=p[0],hi=p[0];
        for (i=1;i<3;++i) {
            if (p[i].x<lo.x) lo.x=p[i].x;
            if (p[i].y<lo.y) lo.y=p[i].y;
            if (p[i].z<lo.z) lo.z=p[i].z;
            if (p[i].x>hi.x) hi.x=p[i].x;
            if (p[i].y>hi.y) hi.y=p[i].y;
            if (p[i].z>hi.z) hi.z=p[i].z;
        }
        return geomlib_aabb_create(lo,hi);
    }
}

static geomlib_aabb merge_boxes(geomlib_aabb a, geomlib_aabb b)
{
    geomlib_vec3 lo,hi;
    lo.x=fminf(a.min.x,b.min.x); lo.y=fminf(a.min.y,b.min.y);
    lo.z=fminf(a.min.z,b.min.z);
    hi.x=fmaxf(a.max.x,b.max.x); hi.y=fmaxf(a.max.y,b.max.y);
    hi.z=fmaxf(a.max.z,b.max.z);
    return geomlib_aabb_create(lo,hi);
}

static int add_node(ov_tree *t, geomlib_aabb box, uint32_t lo, uint32_t hi,
                    int left, int right)
{
    if (t->count==t->cap) {
        size_t nc=t->cap?t->cap*2:64;
        ov_node *nn=(ov_node*)realloc(t->nodes,nc*sizeof(ov_node));
        if (!nn) return -1;
        t->nodes=nn; t->cap=nc;
    }
    t->nodes[t->count].box=box;
    t->nodes[t->count].lo=lo;
    t->nodes[t->count].hi=hi;
    t->nodes[t->count].left=left;
    t->nodes[t->count].right=right;
    return (int)t->count++;
}

/* Split by index range; interior nodes store (split-1, split) as a sentinel
 * so leaves are identified by lo==hi. */
static int build(ov_tree *t, const geomlib_trimesh *m,
                 uint32_t lo, uint32_t hi)
{
    geomlib_aabb box=tri_box(m,lo);
    if (hi==lo) return add_node(t,box,lo,hi,-1,-1);
    {
        uint32_t split=lo+(hi-lo+1)/2;
        int left=build(t,m,lo,split-1);
        int right=build(t,m,split,hi);
        if (left<0||right<0) return -1;
        box=merge_boxes(t->nodes[left].box,t->nodes[right].box);
        return add_node(t,box,split-1,split,left,right);
    }
}

static int append_pair(geomlib_overlap_pairs *p, uint32_t a, uint32_t b)
{
    if (p->count==p->capacity) {
        size_t nc=p->capacity?p->capacity*2:64;
        uint32_t *np=(uint32_t*)realloc(p->pairs,nc*2*sizeof(uint32_t));
        if (!np) return -1;
        p->pairs=np; p->capacity=nc;
    }
    p->pairs[p->count*2]=a;
    p->pairs[p->count*2+1]=b;
    p->count++;
    return 0;
}

static int query(const ov_tree *ta, int ia, const ov_tree *tb, int ib,
                 const geomlib_trimesh *a, const geomlib_trimesh *b,
                 geomlib_overlap_pairs *out)
{
    const ov_node *na,*nb;
    if (ia<0||ib<0) return 0;
    if (!boxes_overlap(ta->nodes[ia].box,tb->nodes[ib].box)) return 0;
    na=&ta->nodes[ia]; nb=&tb->nodes[ib];
    if (na->lo==na->hi && nb->lo==nb->hi)
        return append_pair(out,na->lo,nb->lo);
    if (na->lo==na->hi) {
        geomlib_aabb box=tri_box(a,na->lo);
        uint32_t i;
        for (i=nb->lo;i<=nb->hi;++i)
            if (boxes_overlap(box,tri_box(b,i)) &&
                append_pair(out,na->lo,i)<0) return -1;
        return 0;
    }
    if (nb->lo==nb->hi) {
        geomlib_aabb box=tri_box(b,nb->lo);
        uint32_t i;
        for (i=na->lo;i<=na->hi;++i)
            if (boxes_overlap(tri_box(a,i),box) &&
                append_pair(out,i,nb->lo)<0) return -1;
        return 0;
    }
    if (query(ta,na->left,tb,nb->left,a,b,out)<0) return -1;
    if (query(ta,na->left,tb,nb->right,a,b,out)<0) return -1;
    if (query(ta,na->right,tb,nb->left,a,b,out)<0) return -1;
    if (query(ta,na->right,tb,nb->right,a,b,out)<0) return -1;
    return 0;
}

geomlib_overlap_pairs geomlib_overlap_pairs_create(void)
{
    geomlib_overlap_pairs p;
    memset(&p,0,sizeof(p));
    return p;
}

void geomlib_overlap_pairs_free(geomlib_overlap_pairs *p)
{
    if (!p) return;
    free(p->pairs);
    memset(p,0,sizeof(*p));
}

int geomlib_aabb_overlap_pairs(const geomlib_trimesh *a,
                               const geomlib_trimesh *b,
                               geomlib_overlap_pairs *out)
{
    ov_tree ta,tb;
    int roota,rootb,rc=GEOMLIB_MESH_ERROR_INVALID;
    memset(&ta,0,sizeof(ta));
    memset(&tb,0,sizeof(tb));
    if (!out||!a||!b||!a->indices||!b->indices||
        !a->triangle_count||!b->triangle_count) return rc;
    memset(out,0,sizeof(*out));
    roota=build(&ta,a,0,(uint32_t)(a->triangle_count-1));
    rootb=build(&tb,b,0,(uint32_t)(b->triangle_count-1));
    if (roota<0||rootb<0) goto cleanup;
    rc=query(&ta,roota,&tb,rootb,a,b,out);
    if (rc<0) rc=GEOMLIB_MESH_ERROR_MEMORY;

cleanup:
    free(ta.nodes); free(tb.nodes);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_overlap_pairs_free(out);
    return rc;
}
