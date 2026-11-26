/*
 * geomlib - 2D constrained Delaunay triangulation (see geomlib_cdt2d.h)
 *
 * Clean-room implementation: start from the Bowyer-Watson Delaunay result,
 * recover each segment by flipping crossed edges (inserting a Steiner point
 * when the quadrilateral is not convex), then restore the empty-circumcircle
 * property without touching constrained edges.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_cdt2d.h"
#include "geomlib/geomlib_delaunay2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_cdt2d geomlib_cdt2d_create(void)
{
    geomlib_cdt2d c;
    memset(&c, 0, sizeof(c));
    return c;
}

void geomlib_cdt2d_free(geomlib_cdt2d *c)
{
    if (!c) return;
    free(c->vertices);
    free(c->triangles);
    free(c->constraints);
    memset(c, 0, sizeof(*c));
}

typedef struct { uint32_t v[3]; } tri;

typedef struct {
    geomlib_vec2 *pts;
    size_t np, pcap;
    tri *tris;
    size_t nt, tcap;
    uint32_t *cons;
    size_t nc;
} state;

static double orient(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x-a.x)*((double)c.y-a.y)
         - ((double)b.y-a.y)*((double)c.x-a.x);
}

static int in_circle(state *s, uint32_t ia, uint32_t ib, uint32_t ic, uint32_t ip)
{
    geomlib_vec2 a=s->pts[ia], b=s->pts[ib], c=s->pts[ic], p=s->pts[ip];
    double ax=a.x, ay=a.y, bx=b.x, by=b.y, cx=c.x, cy=c.y, px=p.x, py=p.y;
    double adx=ax-px, ady=ay-py, bdx=bx-px, bdy=by-py, cdx=cx-px, cdy=cy-py;
    double det = (adx*adx+ady*ady)*(bdx*cdy-bdy*cdx)
               - (bdx*bdx+bdy*bdy)*(adx*cdy-ady*cdx)
               + (cdx*cdx+cdy*cdy)*(adx*bdy-ady*bdx);
    double o = orient(a,b,c);
    return o >= 0.0 ? det > 1e-12 : det < -1e-12;
}

static int add_tri(state *s, uint32_t a, uint32_t b, uint32_t c)
{
    if (s->nt == s->tcap) {
        size_t nc = s->tcap ? s->tcap*2 : 64;
        tri *nt = (tri *)realloc(s->tris, nc*sizeof(tri));
        if (!nt) return -1;
        s->tris = nt; s->tcap = nc;
    }
    s->tris[s->nt].v[0]=a; s->tris[s->nt].v[1]=b; s->tris[s->nt].v[2]=c;
    s->nt++;
    return 0;
}

static int add_point(state *s, geomlib_vec2 p)
{
    if (s->np == s->pcap) {
        size_t nc = s->pcap ? s->pcap*2 : 64;
        geomlib_vec2 *np = (geomlib_vec2 *)realloc(s->pts, nc*sizeof(geomlib_vec2));
        if (!np) return -1;
        s->pts = np; s->pcap = nc;
    }
    s->pts[s->np++] = p;
    return (int)(s->np-1);
}

/* Locate edge (a,b): returns the two incident triangle indices and their
 * third vertices. */
static size_t find_edge(state *s, uint32_t a, uint32_t b,
                        size_t *t1, size_t *t2, uint32_t *e1, uint32_t *e2)
{
    size_t n = 0, i;
    *t1 = *t2 = (size_t)-1;
    for (i = 0; i < s->nt; ++i) {
        uint32_t x=s->tris[i].v[0], y=s->tris[i].v[1], z=s->tris[i].v[2];
        uint32_t other;
        if      ((x==a&&y==b)||(x==b&&y==a)) other=z;
        else if ((x==a&&z==b)||(x==b&&z==a)) other=y;
        else if ((y==a&&z==b)||(y==b&&z==a)) other=x;
        else continue;
        if (n == 0) { *t1=i; *e1=other; }
        else { *t2=i; *e2=other; }
        n++;
        if (n == 2) break;
    }
    return n;
}

static int constrained(state *s, uint32_t a, uint32_t b)
{
    size_t i;
    for (i = 0; i+1 < s->nc; i += 2)
        if ((s->cons[i]==a && s->cons[i+1]==b) ||
            (s->cons[i]==b && s->cons[i+1]==a)) return 1;
    return 0;
}

static int proper_cross(geomlib_vec2 a, geomlib_vec2 b,
                        geomlib_vec2 c, geomlib_vec2 d)
{
    double o1=orient(a,b,c), o2=orient(a,b,d);
    double o3=orient(c,d,a), o4=orient(c,d,b);
    return o1*o2 < -1e-14 && o3*o4 < -1e-14;
}

static geomlib_vec2 segment_intersection(geomlib_vec2 a, geomlib_vec2 b,
                                         geomlib_vec2 c, geomlib_vec2 d)
{
    double r0=(double)b.x-a.x, r1=(double)b.y-a.y;
    double s0=(double)d.x-c.x, s1=(double)d.y-c.y;
    double den = r0*s1 - r1*s0;
    double t;
    if (fabs(den) < 1e-18) return geomlib_vec2_create((a.x+b.x)*0.5f,(a.y+b.y)*0.5f);
    t = (((double)c.x-a.x)*s1 - ((double)c.y-a.y)*s0)/den;
    return geomlib_vec2_create((float)((double)a.x+t*r0),
                               (float)((double)a.y+t*r1));
}

/* Flip shared edge using the two incident triangles t1,t2. */
static int flip(state *s, size_t t1, size_t t2, uint32_t c, uint32_t d,
                uint32_t e, uint32_t f)
{
    double oe = orient(s->pts[e], s->pts[f], s->pts[c]);
    double od = orient(s->pts[e], s->pts[f], s->pts[d]);
    if (fabs(oe) < 1e-12 || fabs(od) < 1e-12 || oe*od >= 0) return 0;
    s->tris[t1].v[0]=e; s->tris[t1].v[1]=f; s->tris[t1].v[2]=c;
    s->tris[t2].v[0]=f; s->tris[t2].v[1]=e; s->tris[t2].v[2]=d;
    return 1;
}

static int find_crossing_edge(state *s, uint32_t a, uint32_t b,
                              uint32_t *c, uint32_t *d)
{
    size_t i;
    int k;
    for (i = 0; i < s->nt; ++i)
        for (k = 0; k < 3; ++k) {
            uint32_t x=s->tris[i].v[k], y=s->tris[i].v[(k+1)%3];
            if (x>y) { uint32_t t=x; x=y; y=t; }
            if (x==a && y==b) continue;
            if (x!=a && x!=b && y!=a && y!=b &&
                proper_cross(s->pts[a],s->pts[b],s->pts[x],s->pts[y])) {
                *c=x; *d=y; return 1;
            }
        }
    return 0;
}

static int recover(state *s, uint32_t a, uint32_t b)
{
    int iter;
    for (iter = 0; iter < 400; ++iter) {
        uint32_t c=0, d=0, e=0, f=0, sidx=0;
        size_t t1, t2;
        geomlib_vec2 crossing;
        if (!find_crossing_edge(s, a, b, &c, &d)) return 1;
        if (find_edge(s, c, d, &t1, &t2, &e, &f) != 2) return 0;
        crossing = segment_intersection(s->pts[a], s->pts[b],
                                       s->pts[c], s->pts[d]);
        if (flip(s, t1, t2, c, d, e, f)) continue;
        /* Non-convex quadrilateral: place a Steiner point on cd and split
         * the two incident triangles. */
        sidx = (uint32_t)add_point(s, crossing);
        if (add_tri(s, c, sidx, e) < 0 || add_tri(s, sidx, d, e) < 0 ||
            add_tri(s, d, sidx, f) < 0 || add_tri(s, sidx, c, f) < 0) return 0;
        /* Remove the two original triangles (swap with tail). */
        if (t1 > t2) { size_t tmp=t1; t1=t2; t2=tmp; }
        s->tris[t1] = s->tris[s->nt-2];
        s->tris[t2] = s->tris[s->nt-1];
        s->nt -= 2;
    }
    return 0;
}

static void restore_delaunay(state *s)
{
    int pass;
    for (pass = 0; pass < 200; ++pass) {
        size_t i;
        int k, changed = 0;
        for (i = 0; i < s->nt; ++i)
            for (k = 0; k < 3; ++k) {
                uint32_t a=s->tris[i].v[k], b=s->tris[i].v[(k+1)%3];
                uint32_t c=0, d=0, e=0, f=0;
                size_t t1, t2;
                if (a > b) { uint32_t t=a; a=b; b=t; }
                if (constrained(s,a,b)) continue;
                if (find_edge(s,a,b,&t1,&t2,&e,&f) != 2) continue;
                c=e; d=f;
                if (in_circle(s,a,b,c,d)) {
                    if (flip(s,t1,t2,a,b,c,d)) { changed=1; break; }
                }
            }
        if (!changed) break;
    }
}

int geomlib_cdt2d_build(const geomlib_vec2 *points, size_t count,
                        const uint32_t *segments, size_t segment_count,
                        geomlib_cdt2d *out)
{
    state st;
    geomlib_vec2 *uniq = NULL;
    size_t nu = 0, i;
    geomlib_delaunay2d_result del = geomlib_delaunay2d_result_create();
    geomlib_cdt2d result = geomlib_cdt2d_create();
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !points || count < 3) return rc;
    memset(&st, 0, sizeof(st));
    uniq = (geomlib_vec2 *)malloc((count + 2*segment_count)*sizeof(geomlib_vec2));
    if (!uniq) goto cleanup;
    for (i = 0; i < count; ++i) {
        size_t j;
        int dup = 0;
        for (j = 0; j < nu; ++j)
            if (fabs(uniq[j].x-points[i].x) < 1e-10 &&
                fabs(uniq[j].y-points[i].y) < 1e-10) { dup=1; break; }
        if (!dup) uniq[nu++] = points[i];
    }
    for (i = 0; i+1 < 2*segment_count; i += 2) {
        size_t j;
        int k;
        for (k = 0; k < 2; ++k) {
            uint32_t id = segments[i+k];
            int dup = (id >= count);
            geomlib_vec2 p = dup ? uniq[0] : points[id];
            if (!dup)
                for (j = 0; j < nu; ++j)
                    if (fabs(uniq[j].x-p.x)<1e-10 && fabs(uniq[j].y-p.y)<1e-10) {
                        dup=1; break;
                    }
            if (!dup) uniq[nu++] = p;
        }
    }
    if (nu < 3) goto cleanup;
    if (geomlib_delaunay2d(uniq, nu, &del) != GEOMLIB_MESH_SUCCESS) goto cleanup;
    if (!del.triangle_count) goto cleanup;

    for (i = 0; i < nu; ++i) if (add_point(&st, uniq[i]) < 0) goto cleanup;
    for (i = 0; i < del.triangle_count; ++i)
        if (add_tri(&st, del.triangles[i*3], del.triangles[i*3+1],
                    del.triangles[i*3+2]) < 0) goto cleanup;

    st.cons = (uint32_t *)malloc((segment_count*2+2)*sizeof(uint32_t));
    if (!st.cons) goto cleanup;
    for (i = 0; i+1 < 2*segment_count; i += 2) {
        uint32_t a = segments[i], b = segments[i+1];
        if (a >= nu || b >= nu || a == b) continue;
        if (!recover(&st, a, b)) goto cleanup;
        st.cons[st.nc++] = a; st.cons[st.nc++] = b;
    }
    restore_delaunay(&st);

    result.vertices = st.pts; st.pts = NULL;
    result.vertex_count = st.np;
    result.triangles = (uint32_t *)malloc(st.nt*3*sizeof(uint32_t));
    if (!result.triangles) goto cleanup;
    for (i = 0; i < st.nt; ++i) {
        result.triangles[i*3]   = st.tris[i].v[0];
        result.triangles[i*3+1] = st.tris[i].v[1];
        result.triangles[i*3+2] = st.tris[i].v[2];
    }
    result.triangle_count = st.nt;
    result.constraints = st.cons; st.cons = NULL;
    result.constraint_count = st.nc;
    *out = result;
    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    free(uniq);
    free(st.pts); free(st.tris); free(st.cons);
    geomlib_delaunay2d_result_free(&del);
    if (rc != GEOMLIB_MESH_SUCCESS) geomlib_cdt2d_free(&result);
    return rc;
}

int geomlib_cdt2d_has_edge(const geomlib_cdt2d *c, uint32_t a, uint32_t b)
{
    size_t i;
    if (!c || !c->triangles) return 0;
    for (i = 0; i < c->triangle_count; ++i) {
        uint32_t x=c->triangles[i*3], y=c->triangles[i*3+1], z=c->triangles[i*3+2];
        if ((x==a&&y==b)||(y==a&&x==b)||(x==a&&z==b)||
            (z==a&&x==b)||(y==a&&z==b)||(z==a&&y==b)) return 1;
    }
    return 0;
}
