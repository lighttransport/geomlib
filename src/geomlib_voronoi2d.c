/*
 * geomlib - 2D Voronoi diagram (see geomlib_voronoi2d.h)
 *
 * Clean-room half-plane intersection: a cell is the bounding box clipped by
 * the perpendicular bisector of the site and every other site
 * ( |x - pi|^2 <= |x - pj|^2  <=>  2(pj-pi).x <= |pj|^2 - |pi|^2 ).
 * Clipping uses Sutherland-Hodgman and preserves the CCW winding of the box.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_voronoi2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_voronoi2d geomlib_voronoi2d_create(void)
{
    geomlib_voronoi2d v;
    memset(&v, 0, sizeof(v));
    return v;
}

void geomlib_voronoi2d_free(geomlib_voronoi2d *v)
{
    if (!v) return;
    free(v->sites);
    free(v->cells);
    free(v->offsets);
    memset(v, 0, sizeof(*v));
}

typedef struct { double x, y; } point2;

static int near_equal(point2 a, point2 b)
{
    return fabs(a.x-b.x) < 1e-10 && fabs(a.y-b.y) < 1e-10;
}

/* Clip a CCW polygon to n.x <= c; returns the new vertex count. */
static size_t clip_halfplane(const point2 *in, size_t nin,
                             double nx, double ny, double c,
                             point2 *out)
{
    size_t nout = 0, i;
    const double eps = 1e-12;
    if (!nin) return 0;
    for (i = 0; i < nin; ++i) {
        point2 a = in[i], b = in[(i+1)%nin];
        double fa = nx*a.x + ny*a.y - c;
        double fb = nx*b.x + ny*b.y - c;
        int ain = fa <= eps, bin = fb <= eps;
        if (ain && (nout == 0 || !near_equal(out[nout-1], a))) out[nout++] = a;
        if (ain != bin) {
            double t = fa / (fa - fb);
            point2 q;
            q.x = a.x + t*(b.x-a.x);
            q.y = a.y + t*(b.y-a.y);
            if (nout == 0 || !near_equal(out[nout-1], q)) out[nout++] = q;
        }
    }
    if (nout > 1 && near_equal(out[0], out[nout-1])) nout--;
    return nout;
}

static void init_box(point2 *box, double minx, double miny,
                     double maxx, double maxy)
{
    box[0].x=minx; box[0].y=miny;
    box[1].x=maxx; box[1].y=miny;
    box[2].x=maxx; box[2].y=maxy;
    box[3].x=minx; box[3].y=maxy;
}

/* Clip the box by all bisectors of site i; returns the polygon vertex count
 * and fills `poly`. Scratch buffers are caller-owned with capacity `cap`. */
static size_t build_cell(size_t i, const geomlib_vec2 *sites, size_t nsites,
                         const point2 *box, point2 *bufa, point2 *bufb,
                         size_t cap, point2 *poly)
{
    size_t na = 4, j;
    memcpy(bufa, box, 4*sizeof(point2));
    for (j = 0; j < nsites; ++j) {
        double px, py, qx, qy, nx, ny, c, nb;
        point2 *tmp;
        if (i == j) continue;
        px = sites[i].x; py = sites[i].y;
        qx = sites[j].x; qy = sites[j].y;
        nx = 2.0*(qx-px); ny = 2.0*(qy-py);
        c = qx*qx + qy*qy - px*px - py*py;
        nb = clip_halfplane(bufa, na, nx, ny, c, bufb);
        tmp = bufa; bufa = bufb; bufb = tmp;
        na = nb;
        if (na < 3 || na > cap) break;
    }
    if (na >= 3 && na <= cap) memcpy(poly, bufa, na*sizeof(point2));
    else na = 0;
    return na;
}

int geomlib_compute_voronoi2d(const geomlib_vec2 *points, size_t count,
                              geomlib_vec2 bmin, geomlib_vec2 bmax,
                              geomlib_voronoi2d *out)
{
    geomlib_voronoi2d result = geomlib_voronoi2d_create();
    uint32_t *keep = NULL;
    size_t *counts = NULL, *offsets = NULL, total = 0, i, j, nsites = 0;
    double minx, miny, maxx, maxy;
    point2 box[4], *bufa = NULL, *bufb = NULL, *poly = NULL;
    size_t cap;

    if (!points || !out || count < 1) return GEOMLIB_MESH_ERROR_INVALID;
    minx = (bmin.x < bmax.x) ? bmin.x : bmax.x;
    maxx = (bmin.x < bmax.x) ? bmax.x : bmin.x;
    miny = (bmin.y < bmax.y) ? bmin.y : bmax.y;
    maxy = (bmin.y < bmax.y) ? bmax.y : bmin.y;
    if (!(maxx > minx) || !(maxy > miny)) return GEOMLIB_MESH_ERROR_INVALID;
    init_box(box, minx, miny, maxx, maxy);

    keep = (uint32_t *)malloc(count*sizeof(uint32_t));
    if (!keep) goto fail;
    for (i = 0; i < count; ++i) {
        int dup = 0;
        for (j = 0; j < nsites; ++j) {
            double dx = (double)points[keep[j]].x - points[i].x;
            double dy = (double)points[keep[j]].y - points[i].y;
            if (dx*dx + dy*dy < 1e-18) { dup = 1; break; }
        }
        if (!dup) keep[nsites++] = (uint32_t)i;
    }
    if (!nsites) goto fail;

    result.sites = (geomlib_vec2 *)malloc(nsites*sizeof(geomlib_vec2));
    counts = (size_t *)malloc(nsites*sizeof(size_t));
    offsets = (size_t *)malloc((nsites+1)*sizeof(size_t));
    if (!result.sites || !counts || !offsets) goto fail;
    for (i = 0; i < nsites; ++i) result.sites[i] = points[keep[i]];
    result.offsets = offsets;

    cap = nsites + 8;
    bufa = (point2 *)malloc(cap*sizeof(point2));
    bufb = (point2 *)malloc(cap*sizeof(point2));
    poly = (point2 *)malloc(cap*sizeof(point2));
    if (!bufa || !bufb || !poly) goto fail;

    for (i = 0; i < nsites; ++i) {
        counts[i] = build_cell(i, result.sites, nsites, box, bufa, bufb, cap, poly);
        offsets[i] = total;
        total += counts[i];
    }
    offsets[nsites] = total;

    result.cells = (geomlib_vec2 *)malloc((total ? total : 1)*sizeof(geomlib_vec2));
    if (!result.cells) goto fail;
    for (i = 0; i < nsites; ++i) {
        size_t n = build_cell(i, result.sites, nsites, box, bufa, bufb, cap, poly);
        size_t k;
        if (n > counts[i]) n = counts[i];
        for (k = 0; k < n; ++k)
            result.cells[offsets[i]+k] = geomlib_vec2_create(
                (float)poly[k].x, (float)poly[k].y);
    }

    result.site_count = nsites;
    result.vertex_count = total;
    result.bounds_min = bmin;
    result.bounds_max = bmax;
    free(keep); free(counts); free(bufa); free(bufb); free(poly);
    *out = result;
    return GEOMLIB_MESH_SUCCESS;

fail:
    free(keep); free(counts); free(bufa); free(bufb); free(poly);
    geomlib_voronoi2d_free(&result);
    return GEOMLIB_MESH_ERROR_MEMORY;
}

size_t geomlib_voronoi2d_locate(const geomlib_voronoi2d *v, geomlib_vec2 q)
{
    size_t i, best = 0;
    double best_d = DBL_MAX;
    if (!v || !v->sites || !v->site_count) return (size_t)-1;
    if (q.x < v->bounds_min.x || q.x > v->bounds_max.x ||
        q.y < v->bounds_min.y || q.y > v->bounds_max.y) return (size_t)-1;
    for (i = 0; i < v->site_count; ++i) {
        double dx = (double)q.x - v->sites[i].x;
        double dy = (double)q.y - v->sites[i].y;
        double d = dx*dx + dy*dy;
        if (d < best_d) { best_d = d; best = i; }
    }
    return best;
}

double geomlib_voronoi2d_cell_area(const geomlib_voronoi2d *v, size_t i)
{
    size_t a, b, n;
    double s = 0.0;
    if (!v || !v->offsets || i >= v->site_count) return 0.0;
    a = v->offsets[i]; b = v->offsets[i+1];
    if (b <= a+2) return 0.0;
    for (n = a; n < b; ++n) {
        size_t nx = (n+1 < b) ? n+1 : a;
        s += (double)v->cells[n].x*v->cells[nx].y
           - (double)v->cells[nx].x*v->cells[n].y;
    }
    return 0.5*fabs(s);
}
