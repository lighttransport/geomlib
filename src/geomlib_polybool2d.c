/*
 * geomlib - 2D polygon Boolean operations (see geomlib_polybool2d.h)
 *
 * Clean-room edge-arrangement implementation: split boundaries at crossings,
 * keep sub-edges according to the Boolean predicate, and trace closed loops.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polybool2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_poly_bool_result geomlib_poly_bool_result_create(void)
{
    geomlib_poly_bool_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_poly_bool_result_free(geomlib_poly_bool_result *r)
{
    if (!r) return;
    free(r->loops);
    free(r->offsets);
    memset(r, 0, sizeof(*r));
}

typedef struct { geomlib_vec2 a, b; int from_a; int used; } edge2;

static double dist2(geomlib_vec2 a, geomlib_vec2 b)
{
    double dx=(double)a.x-b.x, dy=(double)a.y-b.y;
    return dx*dx+dy*dy;
}

static double orient(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x-a.x)*((double)c.y-a.y)
         - ((double)b.y-a.y)*((double)c.x-a.x);
}

static double seg_param(geomlib_vec2 p, geomlib_vec2 q, geomlib_vec2 r)
{
    double vx=(double)q.x-p.x, vy=(double)q.y-p.y;
    double l2=vx*vx+vy*vy;
    if (!(l2 > 1e-24)) return 0.0;
    return (((double)r.x-p.x)*vx + ((double)r.y-p.y)*vy)/l2;
}

static double point_line_dist(geomlib_vec2 p, geomlib_vec2 q, geomlib_vec2 r)
{
    double d = fabs(orient(p,q,r));
    double l = hypot((double)q.x-p.x,(double)q.y-p.y);
    return l > 0 ? d/l : d;
}

static int point_in_polygon(geomlib_vec2 q, const geomlib_vec2 *p, size_t n)
{
    size_t i;
    int inside = 0;
    for (i = 0; i < n; ++i) {
        size_t j=(i+1)%n;
        if (((p[i].y>q.y)!=(p[j].y>q.y)) &&
            ((double)q.x < ((double)p[j].x-p[i].x)*((double)q.y-p[i].y)/
                           ((double)p[j].y-p[i].y)+p[i].x))
            inside = !inside;
    }
    return inside;
}

static double polygon_area(const geomlib_vec2 *p, size_t n)
{
    size_t i;
    double s=0.0;
    for (i=0;i<n;++i) {
        size_t j=(i+1)%n;
        s += (double)p[i].x*p[j].y - (double)p[j].x*p[i].y;
    }
    return 0.5*s;
}

static int cmp_double(const void *pa, const void *pb)
{
    double a=*(const double*)pa, b=*(const double*)pb;
    return a<b?-1:(a>b?1:0);
}

static geomlib_vec2 lerp(geomlib_vec2 a, geomlib_vec2 b, double t)
{
    return geomlib_vec2_create((float)(a.x+t*(b.x-a.x)),
                               (float)(a.y+t*(b.y-a.y)));
}

static int add_edge(edge2 **arr, size_t *n, size_t *cap,
                    geomlib_vec2 a, geomlib_vec2 b, int from_a, double eps2)
{
    double dx=(double)b.x-a.x, dy=(double)b.y-a.y;
    if (dx*dx+dy*dy < eps2) return 0;
    if (*n==*cap) {
        size_t nc=*cap?(*cap)*2:64;
        edge2 *na=(edge2*)realloc(*arr,nc*sizeof(edge2));
        if (!na) return -1;
        *arr=na; *cap=nc;
    }
    (*arr)[*n].a=a; (*arr)[*n].b=b;
    (*arr)[*n].from_a=from_a; (*arr)[*n].used=0;
    (*n)++;
    return 0;
}

/* Split edges of `src` at intersections with edges of `other`. */
static int split_edges(const geomlib_vec2 *src, size_t ns,
                       const geomlib_vec2 *other, size_t no,
                       int from_a, edge2 **edges, size_t *ne,
                       size_t *cap, double eps)
{
    double *ts;
    size_t i, j;
    ts=(double*)malloc((no+2)*sizeof(double));
    if (!ts) return -1;
    for (i=0;i<ns;++i) {
        size_t j2=(i+1)%ns, nt=0;
        geomlib_vec2 p=src[i], q=src[j2];
        ts[nt++]=0.0; ts[nt++]=1.0;
        for (j=0;j<no;++j) {
            size_t k=(j+1)%no;
            geomlib_vec2 r=other[j], s=other[k];
            double d1=orient(p,q,r), d2=orient(p,q,s);
            double d3=orient(r,s,p), d4=orient(r,s,q);
            if (d1*d2<0.0 && d3*d4<0.0) {
                double t = d3/(d3-d4);
                if (t>eps && t<1.0-eps) ts[nt++]=t;
            }
            if (point_line_dist(r,p,q)<eps) {
                double t=seg_param(p,q,r);
                if (t>eps && t<1.0-eps) ts[nt++]=t;
            }
            if (point_line_dist(s,p,q)<eps) {
                double t=seg_param(p,q,s);
                if (t>eps && t<1.0-eps) ts[nt++]=t;
            }
        }
        qsort(ts,nt,sizeof(double),cmp_double);
        for (j=0;j+1<nt;++j) {
            if (ts[j+1]-ts[j] < eps) continue;
            if (add_edge(edges,ne,cap,lerp(p,q,ts[j]),lerp(p,q,ts[j+1]),
                         from_a, eps*eps)<0) { free(ts); return -1; }
        }
    }
    free(ts);
    return 0;
}

int geomlib_polygon_boolean(const geomlib_vec2 *a, size_t na,
                            const geomlib_vec2 *b, size_t nb,
                            geomlib_poly_bool_op op,
                            geomlib_poly_bool_result *out)
{
    geomlib_poly_bool_result result = geomlib_poly_bool_result_create();
    geomlib_vec2 *pa=NULL, *pb=NULL;
    edge2 *edges=NULL;
    size_t ne=0, cap=0, i;
    double eps, eps2, extent=0.0;
    geomlib_vec2 lo, hi;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !a || !b || na<3 || nb<3) return rc;
    lo=a[0]; hi=a[0];
    for (i=0;i<na;++i) {
        if (a[i].x<lo.x) lo.x=a[i].x;
        if (a[i].y<lo.y) lo.y=a[i].y;
        if (a[i].x>hi.x) hi.x=a[i].x;
        if (a[i].y>hi.y) hi.y=a[i].y;
    }
    for (i=0;i<nb;++i) {
        if (b[i].x<lo.x) lo.x=b[i].x;
        if (b[i].y<lo.y) lo.y=b[i].y;
        if (b[i].x>hi.x) hi.x=b[i].x;
        if (b[i].y>hi.y) hi.y=b[i].y;
    }
    extent = fmax((double)hi.x-lo.x,(double)hi.y-lo.y);
    if (!(extent > 1e-9)) return rc;
    eps = extent*1e-7;
    eps2 = eps*eps;

    pa=(geomlib_vec2*)malloc(na*sizeof(geomlib_vec2));
    pb=(geomlib_vec2*)malloc(nb*sizeof(geomlib_vec2));
    if (!pa||!pb) goto cleanup;
    /* Normalize to CCW so the keep rules have consistent orientation. */
    if (polygon_area(a,na) < 0.0)
        for (i=0;i<na;++i) pa[i]=a[na-1-i];
    else memcpy(pa,a,na*sizeof(geomlib_vec2));
    if (polygon_area(b,nb) < 0.0)
        for (i=0;i<nb;++i) pb[i]=b[nb-1-i];
    else memcpy(pb,b,nb*sizeof(geomlib_vec2));

    if (split_edges(pa,na,pb,nb,1,&edges,&ne,&cap,eps)<0 ||
        split_edges(pb,nb,pa,na,0,&edges,&ne,&cap,eps)<0) goto cleanup;

    /* Drop sub-edges not selected by the operation. */
    {
        size_t w=0;
        for (i=0;i<ne;++i) {
            geomlib_vec2 mid=lerp(edges[i].a,edges[i].b,0.5);
            int in_a=point_in_polygon(mid,pa,na);
            int in_b=point_in_polygon(mid,pb,nb);
            int keep;
            if (edges[i].from_a) keep = op==GEOMLIB_POLY_UNION ? !in_b
                                 : op==GEOMLIB_POLY_INTERSECTION ? in_b : !in_b;
            else keep = op==GEOMLIB_POLY_UNION ? !in_a
                      : op==GEOMLIB_POLY_INTERSECTION ? in_a : in_a;
            if (!keep) continue;
            if (op==GEOMLIB_POLY_DIFFERENCE && !edges[i].from_a) {
                geomlib_vec2 t=edges[i].a;
                edges[i].a=edges[i].b; edges[i].b=t;
            }
            if (w!=i) edges[w]=edges[i];
            w++;
        }
        ne=w;
    }

    /* Trace loops by matching edge endpoints. */
    {
        geomlib_vec2 *verts=NULL;
        size_t *offsets=NULL, nv=0, nloops=0, ocap=0;
        offsets=(size_t*)malloc((ne+1)*sizeof(size_t));
        if (!offsets) goto cleanup;
        verts=(geomlib_vec2*)malloc((ne*2+2)*sizeof(geomlib_vec2));
        if (!verts) { free(offsets); goto cleanup; }
        for (i=0;i<ne;++i) {
            geomlib_vec2 start, cur;
            size_t count=0;
            if (edges[i].used) continue;
            start=cur=edges[i].a;
            for (;;) {
                size_t j;
                int found=0;
                verts[nv++]=cur; count++;
                edges[i].used=1;
                cur=edges[i].b;
                if (dist2(cur,start) < eps2) break;
                for (j=0;j<ne;++j) {
                    if (edges[j].used) continue;
                    if (dist2(edges[j].a,cur)<eps2) {
                        i=j; found=1; break;
                    }
                }
                if (!found) break;
            }
            if (count>=3) {
                if (nloops==ocap) {
                    size_t nc=ocap?ocap*2:8;
                    size_t *no=(size_t*)realloc(offsets,nc*sizeof(size_t));
                    if (!no) { free(verts); free(offsets); goto cleanup; }
                    offsets=no; ocap=nc;
                }
                offsets[nloops++]=nv-count;
            } else nv-=count;
        }
        offsets[nloops]=nv;
        result.loops=verts;
        result.offsets=offsets;
        result.loop_count=nloops;
    }

    *out=result;
    rc=GEOMLIB_MESH_SUCCESS;
cleanup:
    free(pa); free(pb); free(edges);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_poly_bool_result_free(&result);
    return rc;
}
