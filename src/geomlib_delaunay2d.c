/*
 * geomlib - 2D Delaunay triangulation (see geomlib_delaunay2d.h)
 *
 * Clean-room Bowyer-Watson incremental insertion. C11, no dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_delaunay2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_delaunay2d_result geomlib_delaunay2d_result_create(void)
{
    geomlib_delaunay2d_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_delaunay2d_result_free(geomlib_delaunay2d_result *r)
{
    if (!r) return;
    free(r->triangles);
    memset(r, 0, sizeof(*r));
}

typedef struct { uint32_t a, b, c; } tri;
typedef struct { uint32_t a, b; } edge;

static int cmp_edge(const void *pa, const void *pb)
{
    const edge *x=(const edge*)pa, *y=(const edge*)pb;
    if (x->a<y->a) return -1;
    if (x->a>y->a) return 1;
    if (x->b<y->b) return -1;
    if (x->b>y->b) return 1;
    return 0;
}

static int in_circumcircle(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c,
                           geomlib_vec2 p)
{
    double ax=a.x, ay=a.y, bx=b.x, by=b.y, cx=c.x, cy=c.y, px=p.x, py=p.y;
    double adx=ax-px, ady=ay-py, bdx=bx-px, bdy=by-py, cdx=cx-px, cdy=cy-py;
    double det = (adx*adx+ady*ady)*(bdx*cdy-bdy*cdx)
               - (bdx*bdx+bdy*bdy)*(adx*cdy-ady*cdx)
               + (cdx*cdx+cdy*cdy)*(adx*bdy-ady*bdx);
    /* Orient so a CCW triangle gives positive inside; the caller may pass any
     * winding, so use the triangle's own orientation sign. */
    double orient = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
    return orient >= 0.0 ? det > 1e-12 : det < -1e-12;
}

int geomlib_delaunay2d(const geomlib_vec2 *points, size_t count,
                       geomlib_delaunay2d_result *out)
{
    geomlib_vec2 *pts;
    double minx=DBL_MAX,miny=DBL_MAX,maxx=-DBL_MAX,maxy=-DBL_MAX;
    double dx,dy,midx,midy;
    tri *tris;
    size_t ntri, cap, i, j, k;
    uint32_t n;

    if (!points || !out || count < 3) return GEOMLIB_MESH_ERROR_INVALID;
    pts = (geomlib_vec2 *)malloc((count+3)*sizeof(geomlib_vec2));
    if (!pts) return GEOMLIB_MESH_ERROR_MEMORY;
    n = 0;
    for (i = 0; i < count; ++i) {
        int dup = 0;
        for (j = 0; j < n; ++j)
            if (fabs(pts[j].x-points[i].x) < 1e-12 &&
                fabs(pts[j].y-points[i].y) < 1e-12) { dup = 1; break; }
        if (dup) continue;
        pts[n] = points[i];
        if (pts[n].x<minx) minx=pts[n].x;
        if (pts[n].y<miny) miny=pts[n].y;
        if (pts[n].x>maxx) maxx=pts[n].x;
        if (pts[n].y>maxy) maxy=pts[n].y;
        n++;
    }
    if (n < 3) { free(pts); return GEOMLIB_MESH_ERROR_INVALID; }
    dx = maxx-minx; dy = maxy-miny;
    if (dx < 1e-12 && dy < 1e-12) { free(pts); return GEOMLIB_MESH_ERROR_INVALID; }
    midx = (minx+maxx)*0.5;
    midy = (miny+maxy)*0.5;
    {
        double d = (dx>dy?dx:dy)*20.0;
        pts[n] = geomlib_vec2_create((float)(midx-2*d), (float)(midy-d));
        pts[n+1] = geomlib_vec2_create((float)midx, (float)(midy+2*d));
        pts[n+2] = geomlib_vec2_create((float)(midx+2*d), (float)(midy-d));
    }

    cap = 8*n+64;
    tris = (tri *)malloc(cap*sizeof(tri));
    if (!tris) { free(pts); return GEOMLIB_MESH_ERROR_MEMORY; }
    ntri = 1;
    tris[0].a=n; tris[0].b=n+1; tris[0].c=n+2;

    for (i = 0; i < n; ++i) {
        uint8_t *bad = (uint8_t *)calloc(ntri, sizeof(uint8_t));
        edge *boundary;
        size_t nbad = 0, nboundary = 0;
        if (!bad) { free(tris); free(pts); return GEOMLIB_MESH_ERROR_MEMORY; }
        for (j = 0; j < ntri; ++j) {
            if (!in_circumcircle(pts[tris[j].a], pts[tris[j].b], pts[tris[j].c], pts[i]))
                continue;
            bad[j] = 1;
            nbad++;
        }
        if (!nbad) { free(bad); continue; }
        boundary = (edge *)malloc(3*nbad*sizeof(edge));
        if (!boundary) { free(bad); free(tris); free(pts); return GEOMLIB_MESH_ERROR_MEMORY; }
        for (j = 0; j < ntri; ++j) {
            if (!bad[j]) continue;
            uint32_t v[3] = {tris[j].a, tris[j].b, tris[j].c};
            for (k = 0; k < 3; ++k) {
                uint32_t a=v[k], b=v[(k+1)%3];
                int shared = 0;
                size_t q;
                for (q = 0; q < ntri; ++q) {
                    if (!bad[q] || q==j) continue;
                    uint32_t w[3] = {tris[q].a, tris[q].b, tris[q].c};
                    int m;
                    for (m = 0; m < 3; ++m)
                        if (((w[m]==a && w[(m+1)%3]==b) ||
                             (w[m]==b && w[(m+1)%3]==a))) { shared=1; break; }
                    if (shared) break;
                }
                if (!shared) {
                    boundary[nboundary].a = a<b?a:b;
                    boundary[nboundary].b = a<b?b:a;
                    nboundary++;
                }
            }
        }
        qsort(boundary, nboundary, sizeof(edge), cmp_edge);
        /* Remove bad triangles (swap with tail). */
        {
            size_t w = 0;
            for (j = 0; j < ntri; ++j) {
                if (bad[j]) continue;
                tris[w++] = tris[j];
            }
            ntri = w;
        }
        free(bad);
        for (j = 0; j < nboundary; ++j) {
            if (ntri >= cap) {
                tri *nt;
                cap *= 2;
                nt = (tri *)realloc(tris, cap*sizeof(tri));
                if (!nt) { free(boundary); free(tris); free(pts); return GEOMLIB_MESH_ERROR_MEMORY; }
                tris = nt;
            }
            tris[ntri].a = (uint32_t)i;
            tris[ntri].b = boundary[j].a;
            tris[ntri].c = boundary[j].b;
            ntri++;
        }
        free(boundary);
    }

    {
        geomlib_delaunay2d_result result = geomlib_delaunay2d_result_create();
        size_t w = 0;
        for (j = 0; j < ntri; ++j) {
            if (tris[j].a>=n || tris[j].b>=n || tris[j].c>=n) continue;
            result.triangles = (uint32_t *)realloc(result.triangles,
                                                  (w+1)*3*sizeof(uint32_t));
            if (!result.triangles) {
                geomlib_delaunay2d_result_free(&result);
                free(tris); free(pts);
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
            result.triangles[w*3]=tris[j].a;
            result.triangles[w*3+1]=tris[j].b;
            result.triangles[w*3+2]=tris[j].c;
            w++;
        }
        result.triangle_count = w;
        *out = result;
    }
    free(tris);
    free(pts);
    return GEOMLIB_MESH_SUCCESS;
}
