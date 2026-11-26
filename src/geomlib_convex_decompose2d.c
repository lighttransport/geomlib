/*
 * geomlib - 2D convex decomposition (see geomlib_convex_decompose2d.h)
 *
 * Clean-room greedy reflex-vertex diagonal cutting. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_convex_decompose2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_poly_decomposition geomlib_poly_decomposition_create(void)
{
    geomlib_poly_decomposition d;
    memset(&d, 0, sizeof(d));
    return d;
}

void geomlib_poly_decomposition_free(geomlib_poly_decomposition *d)
{
    if (!d) return;
    free(d->parts);
    free(d->offsets);
    memset(d, 0, sizeof(*d));
}

static double cross(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x-a.x)*((double)c.y-a.y)
         - ((double)b.y-a.y)*((double)c.x-a.x);
}

static int is_reflex(const geomlib_vec2 *p, size_t n, size_t i)
{
    size_t a=(i+n-1)%n, b=(i+1)%n;
    return cross(p[a],p[i],p[b]) < -1e-10;
}

static int all_convex(const geomlib_vec2 *p, size_t n)
{
    size_t i;
    for (i=0;i<n;++i) if (is_reflex(p,n,i)) return 0;
    return 1;
}

static int proper_cross(geomlib_vec2 a, geomlib_vec2 b,
                        geomlib_vec2 c, geomlib_vec2 d)
{
    double o1=cross(a,b,c),o2=cross(a,b,d);
    double o3=cross(c,d,a),o4=cross(c,d,b);
    return o1*o2 < -1e-12 && o3*o4 < -1e-12;
}

/* Is segment p[i]-p[j] a valid internal diagonal? */
static int valid_diagonal(const geomlib_vec2 *p, size_t n, size_t i, size_t j)
{
    size_t k;
    geomlib_vec2 mid=geomlib_vec2_create((p[i].x+p[j].x)*0.5f,
                                         (p[i].y+p[j].y)*0.5f);
    if (j==(i+1)%n || i==(j+1)%n) return 0;
    for (k=0;k<n;++k) {
        size_t l=(k+1)%n;
        if (k==i||k==j||l==i||l==j) continue;
        if (proper_cross(p[i],p[j],p[k],p[l])) return 0;
    }
    /* Midpoint must be inside (ray casting). */
    {
        int inside=0;
        for (k=0;k<n;++k) {
            size_t l=(k+1)%n;
            if (((p[k].y>mid.y)!=(p[l].y>mid.y)) &&
                ((double)mid.x < ((double)p[l].x-p[k].x)*((double)mid.y-p[k].y)/
                                ((double)p[l].y-p[k].y)+p[k].x))
                inside=!inside;
        }
        if (!inside) return 0;
    }
    return 1;
}

static size_t count_reflex(const geomlib_vec2 *p, size_t n)
{
    size_t i,c=0;
    for (i=0;i<n;++i) if (is_reflex(p,n,i)) c++;
    return c;
}

typedef struct {
    geomlib_vec2 *parts;
    size_t *offsets;
    size_t nparts, cap, ocap;
} builder;

static int add_part(builder *b, const geomlib_vec2 *p, size_t n)
{
    size_t i;
    if (b->nparts==b->ocap) {
        size_t nc=b->ocap?b->ocap*2:16;
        size_t *no=(size_t*)realloc(b->offsets,nc*sizeof(size_t));
        if (!no) return -1;
        b->offsets=no; b->ocap=nc;
    }
    if (b->offsets[b->nparts]+n > b->cap) {
        size_t need=b->offsets[b->nparts]+n;
        size_t nc=b->cap?b->cap*2:64;
        geomlib_vec2 *np;
        while (nc<need) nc*=2;
        np=(geomlib_vec2*)realloc(b->parts,nc*sizeof(geomlib_vec2));
        if (!np) return -1;
        b->parts=np; b->cap=nc;
    }
    for (i=0;i<n;++i) b->parts[b->offsets[b->nparts]+i]=p[i];
    b->nparts++;
    b->offsets[b->nparts]=b->offsets[b->nparts-1]+n;
    return 0;
}

static int decompose_rec(builder *b, const geomlib_vec2 *p, size_t n, int depth)
{
    size_t i,j,bi=(size_t)-1,bj=(size_t)-1;
    size_t best=count_reflex(p,n);
    geomlib_vec2 *a,*c;
    if (n<3) return GEOMLIB_MESH_ERROR_INVALID;
    if (all_convex(p,n) || depth>32) return add_part(b,p,n);
    for (i=0;i<n;++i) {
        if (!is_reflex(p,n,i)) continue;
        for (j=0;j<n;++j) {
            size_t na,nb,ra,rb;
            if (!valid_diagonal(p,n,i,j)) continue;
            na=(j>=i)?(j-i+1):(n-i+j+1);
            nb=n+2-na;
            if (na<3||nb<3||na>256||nb>256) continue;
            ra=rb=0;
            {
                geomlib_vec2 tmp[256];
                size_t k,q=0;
                for (k=i;;k=(k+1)%n) {
                    tmp[q++]=p[k];
                    if (k==j) break;
                }
                ra=count_reflex(tmp,q);
                q=0;
                for (k=j;;k=(k+1)%n) {
                    tmp[q++]=p[k];
                    if (k==i) break;
                }
                rb=count_reflex(tmp,q);
            }
            if (ra+rb<best) { best=ra+rb; bi=i; bj=j; }
        }
    }
    if (bi==(size_t)-1) return add_part(b,p,n);
    a=(geomlib_vec2*)malloc(256*sizeof(geomlib_vec2));
    c=(geomlib_vec2*)malloc(256*sizeof(geomlib_vec2));
    if (!a||!c) { free(a); free(c); return GEOMLIB_MESH_ERROR_MEMORY; }
    {
        size_t na=0,nb=0,k;
        for (k=bi;;k=(k+1)%n) { a[na++]=p[k]; if (k==bj) break; }
        for (k=bj;;k=(k+1)%n) { c[nb++]=p[k]; if (k==bi) break; }
        {
            int rc=decompose_rec(b,a,na,depth+1);
            if (rc==GEOMLIB_MESH_SUCCESS) rc=decompose_rec(b,c,nb,depth+1);
            free(a); free(c);
            return rc;
        }
    }
}

int geomlib_convex_decompose2d(const geomlib_vec2 *polygon, size_t n,
                               geomlib_poly_decomposition *out)
{
    geomlib_poly_decomposition result=geomlib_poly_decomposition_create();
    builder b;
    int rc;
    memset(&b,0,sizeof(b));
    b.ocap=4;
    b.offsets=(size_t*)malloc(b.ocap*sizeof(size_t));
    if (!b.offsets) return GEOMLIB_MESH_ERROR_MEMORY;
    memset(b.offsets,0,b.ocap*sizeof(size_t));
    if (!out || !polygon || n<3) { free(b.offsets); return GEOMLIB_MESH_ERROR_INVALID; }
    rc=decompose_rec(&b,polygon,n,0);
    if (rc==GEOMLIB_MESH_SUCCESS) {
        result.parts=b.parts; b.parts=NULL;
        result.offsets=b.offsets; b.offsets=NULL;
        result.part_count=b.nparts;
        *out=result;
    }
    free(b.parts); free(b.offsets);
    return rc;
}
