/*
 * geomlib - 2D medial axis (see geomlib_medial2d.h)
 *
 * Clean-room sampled-boundary Voronoi construction: Voronoi vertices of
 * boundary samples are candidate centres of maximal inscribed circles.
 * Interior candidates supported by at least two different boundary edges are
 * kept and joined along Voronoi edges.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_medial2d.h"
#include "geomlib/geomlib_voronoi2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_medial2d geomlib_medial2d_create(void)
{
    geomlib_medial2d m;
    memset(&m, 0, sizeof(m));
    return m;
}

void geomlib_medial2d_free(geomlib_medial2d *m)
{
    if (!m) return;
    free(m->segments);
    memset(m, 0, sizeof(*m));
}

typedef struct { geomlib_vec2 p; uint32_t e1, e2; } sample;

typedef struct {
    geomlib_vec2 p;
    uint32_t sites[16];
    size_t nsites;
    int medial;
} corner;

typedef struct { size_t a, b; } link;

static double dist2(geomlib_vec2 a, geomlib_vec2 b)
{
    double dx=(double)a.x-b.x, dy=(double)a.y-b.y;
    return dx*dx+dy*dy;
}

static int point_in_polygon(geomlib_vec2 q, const geomlib_vec2 *p, size_t n)
{
    size_t i;
    int inside = 0;
    for (i = 0; i < n; ++i) {
        size_t j = (i+1)%n;
        if (((p[i].y > q.y) != (p[j].y > q.y)) &&
            ((double)q.x < ((double)p[j].x-p[i].x)*((double)q.y-p[i].y) /
                           ((double)p[j].y-p[i].y) + p[i].x))
            inside = !inside;
    }
    return inside;
}

static double point_seg_dist2(geomlib_vec2 q, geomlib_vec2 a, geomlib_vec2 b)
{
    double vx=(double)b.x-a.x, vy=(double)b.y-a.y;
    double wx=(double)q.x-a.x, wy=(double)q.y-a.y;
    double len2 = vx*vx+vy*vy, t;
    geomlib_vec2 c;
    if (len2 < 1e-20) return dist2(q,a);
    t = (wx*vx+wy*vy)/len2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    c = geomlib_vec2_create((float)(a.x+t*(b.x-a.x)), (float)(a.y+t*(b.y-a.y)));
    return dist2(q,c);
}

static int strictly_inside(geomlib_vec2 q, const geomlib_vec2 *p, size_t n,
                           double eps2)
{
    size_t i;
    if (!point_in_polygon(q,p,n)) return 0;
    for (i = 0; i < n; ++i)
        if (point_seg_dist2(q,p[i],p[(i+1)%n]) < eps2) return 0;
    return 1;
}

static size_t find_or_add_corner(corner **cs, size_t *n, size_t *cap,
                                 geomlib_vec2 p, double eps2)
{
    size_t i;
    for (i = 0; i < *n; ++i)
        if (dist2((*cs)[i].p,p) < eps2) return i;
    if (*n == *cap) {
        size_t nc = *cap ? (*cap)*2 : 64;
        corner *ncorners = (corner *)realloc(*cs, nc*sizeof(corner));
        if (!ncorners) return (size_t)-1;
        *cs = ncorners; *cap = nc;
    }
    memset(&(*cs)[*n], 0, sizeof(corner));
    (*cs)[*n].p = p;
    (*n)++;
    return *n-1;
}

static void add_site(corner *c, uint32_t s)
{
    size_t i;
    for (i = 0; i < c->nsites; ++i) if (c->sites[i] == s) return;
    if (c->nsites < 16) c->sites[c->nsites++] = s;
}

static int add_link(link **ls, size_t *n, size_t *cap, size_t a, size_t b)
{
    size_t i;
    if (a == b) return 0;
    for (i = 0; i < *n; ++i)
        if (((*ls)[i].a==a&&(*ls)[i].b==b)||((*ls)[i].a==b&&(*ls)[i].b==a))
            return 0;
    if (*n == *cap) {
        size_t nc = *cap ? (*cap)*2 : 128;
        link *nlinks = (link *)realloc(*ls, nc*sizeof(link));
        if (!nlinks) return -1;
        *ls = nlinks; *cap = nc;
    }
    (*ls)[*n].a=a; (*ls)[*n].b=b; (*n)++;
    return 0;
}

int geomlib_medial_axis2d(const geomlib_vec2 *polygon, size_t count,
                          float spacing, geomlib_medial2d *out)
{
    geomlib_medial2d result = geomlib_medial2d_create();
    sample *samples = NULL;
    geomlib_vec2 *pts = NULL;
    corner *corners = NULL;
    link *links = NULL;
    geomlib_voronoi2d vor = geomlib_voronoi2d_create();
    size_t ns=0, scap=0, nc=0, ccap=0, nl=0, lcap=0, i, k;
    double minx=DBL_MAX,miny=DBL_MAX,maxx=-DBL_MAX,maxy=-DBL_MAX,diag,eps2;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !polygon || count < 3) return rc;
    for (i = 0; i < count; ++i) {
        if (polygon[i].x<minx) minx=polygon[i].x;
        if (polygon[i].y<miny) miny=polygon[i].y;
        if (polygon[i].x>maxx) maxx=polygon[i].x;
        if (polygon[i].y>maxy) maxy=polygon[i].y;
    }
    diag = hypot(maxx-minx, maxy-miny);
    if (!(diag > 1e-9)) return rc;
    if (!(spacing > 0.0f)) spacing = (float)(0.02*diag);
    eps2 = diag*diag*1e-10;

    for (k = 0; k < count; ++k) {
        geomlib_vec2 a=polygon[k], b=polygon[(k+1)%count];
        double len = hypot((double)b.x-a.x,(double)b.y-a.y);
        int parts = (int)ceil(len/spacing);
        int t;
        if (parts < 1) parts = 1;
        for (t = 0; t < parts; ++t) {
            double f = (double)t/parts;
            if (ns == scap) {
                size_t ncap = scap ? scap*2 : 128;
                sample *nsamples = (sample *)realloc(samples, ncap*sizeof(sample));
                geomlib_vec2 *npts;
                if (!nsamples) goto cleanup;
                samples = nsamples;
                npts = (geomlib_vec2 *)realloc(pts, ncap*sizeof(geomlib_vec2));
                if (!npts) goto cleanup;
                pts = npts;
                scap = ncap;
            }
            samples[ns].p = geomlib_vec2_create((float)(a.x+f*(b.x-a.x)),
                                                (float)(a.y+f*(b.y-a.y)));
            pts[ns] = samples[ns].p;
            samples[ns].e1 = (uint32_t)k;
            samples[ns].e2 = t==0 ? (uint32_t)((k+count-1)%count) : UINT32_MAX;
            ns++;
        }
    }
    if (ns < 4) goto cleanup;

    {
        geomlib_vec2 lo = geomlib_vec2_create((float)minx,(float)miny);
        geomlib_vec2 hi = geomlib_vec2_create((float)maxx,(float)maxy);
        if (geomlib_compute_voronoi2d(pts, ns, lo, hi, &vor) !=
            GEOMLIB_MESH_SUCCESS) goto cleanup;
    }

    for (i = 0; i < vor.site_count; ++i) {
        size_t a=vor.offsets[i], b=vor.offsets[i+1], j;
        for (j = a; j < b; ++j) {
            size_t na = (j+1 < b) ? j+1 : a;
            size_t cj = find_or_add_corner(&corners,&nc,&ccap,vor.cells[j],eps2);
            size_t cn = find_or_add_corner(&corners,&nc,&ccap,vor.cells[na],eps2);
            if (cj==(size_t)-1 || cn==(size_t)-1) goto cleanup;
            add_site(&corners[cj], (uint32_t)i);
            add_site(&corners[cn], (uint32_t)i);
            if (add_link(&links,&nl,&lcap,cj,cn) < 0) goto cleanup;
        }
    }

    for (i = 0; i < nc; ++i) {
        corner *c = &corners[i];
        size_t edges16[16];
        size_t ne = 0, j, q;
        if (!strictly_inside(c->p, polygon, count, eps2)) continue;
        for (j = 0; j < c->nsites; ++j) {
            uint32_t s = c->sites[j];
            uint32_t e[2] = {samples[s].e1, samples[s].e2};
            int t;
            for (t = 0; t < 2; ++t) {
                if (e[t] == UINT32_MAX) continue;
                for (q = 0; q < ne; ++q) if (edges16[q]==e[t]) break;
                if (q == ne && ne < 16) edges16[ne++] = e[t];
            }
        }
        if (ne >= 2) c->medial = 1;
    }

    result.segments = (geomlib_vec2 *)malloc((nl ? nl : 1)*2*sizeof(geomlib_vec2));
    if (!result.segments) goto cleanup;
    for (i = 0; i < nl; ++i) {
        if (!corners[links[i].a].medial || !corners[links[i].b].medial) continue;
        result.segments[result.segment_count*2] = corners[links[i].a].p;
        result.segments[result.segment_count*2+1] = corners[links[i].b].p;
        result.segment_count++;
    }

    *out = result;
    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    free(samples); free(pts); free(corners); free(links);
    geomlib_voronoi2d_free(&vor);
    if (rc != GEOMLIB_MESH_SUCCESS) geomlib_medial2d_free(&result);
    return rc;
}
