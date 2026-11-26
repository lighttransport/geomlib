/*
 * geomlib - MMP continuous geodesic (see geomlib_geodesic_mmp.h)
 *
 * Clean-room window propagation: each window stores an interval on an edge and
 * the unfolded source position in the canonical plane of a face incident to
 * that edge. Windows are expanded across the opposite face in increasing
 * distance order. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_mmp.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MMP_MAX_WINDOWS 40000
#define MMP_EPS 1e-10
#define MMP_MAX_DEPTH 64

typedef struct {
    uint32_t u, v;      /* directed edge carrying the window */
    uint32_t face;      /* face whose canonical frame stores `sx,sy` */
    double t0, t1;      /* interval along the edge */
    double sx, sy;      /* unfolded source in that frame */
    double key;         /* minimum distance over the interval */
    int processed;
    int depth;
} mmp_window;

typedef struct {
    uint32_t a, b;      /* undirected edge key (a < b) */
    uint32_t f0, f1;    /* incident faces, f1 == UINT32_MAX when absent */
} mmp_edge;

static double edge_len(const geomlib_trimesh *m, uint32_t a, uint32_t b)
{
    geomlib_vec3 d = geomlib_vec3_sub(m->vertices[a], m->vertices[b]);
    return sqrt((double)geomlib_vec3_dot(d, d));
}

static void third_vertex(const geomlib_trimesh *m, uint32_t f,
                         uint32_t u, uint32_t v, uint32_t *w)
{
    uint32_t a = m->indices[f * 3], b = m->indices[f * 3 + 1],
             c = m->indices[f * 3 + 2];
    *w = (a != u && a != v) ? a : ((b != u && b != v) ? b : c);
}

static int find_edge(const mmp_edge *edges, size_t ne, uint32_t a, uint32_t b)
{
    size_t i;
    if (a > b) { uint32_t t = a; a = b; b = t; }
    for (i = 0; i < ne; ++i)
        if (edges[i].a == a && edges[i].b == b) return (int)i;
    return -1;
}

static double window_key(const mmp_window *w)
{
    double px = w->sx;
    double t = px < w->t0 ? w->t0 : (px > w->t1 ? w->t1 : px);
    double dx = t - px;
    return sqrt(dx * dx + w->sy * w->sy);
}

static double dist_at(const mmp_window *w, double t)
{
    double dx = t - w->sx;
    return sqrt(dx * dx + w->sy * w->sy);
}

static int add_window(mmp_window **arr, size_t *n, size_t *cap, mmp_window w)
{
    if (*n == *cap) {
        size_t nc = *cap ? *cap * 2 : 128;
        mmp_window *p = (mmp_window *)realloc(*arr, nc * sizeof(mmp_window));
        if (!p) return 0;
        *arr = p;
        *cap = nc;
    }
    w.key = window_key(&w);
    w.processed = 0;
    if (w.depth > MMP_MAX_DEPTH) return 1;
    (*arr)[(*n)++] = w;
    return 1;
}

/* Rotate a 2D point by `ang`. */
static void rot(double x, double y, double ang, double *ox, double *oy)
{
    double c = cos(ang), s = sin(ang);
    *ox = c * x - s * y;
    *oy = s * x + c * y;
}

/* Propagate a window from edge (u,v) of face `face` across one of the two
 * other edges to the third vertex `w`. The source is already expressed in
 * `face`'s canonical frame for edge (u,v); after propagation the new window is
 * expressed in that same face's canonical frame for the new edge. */
static int propagate(const geomlib_trimesh *m, mmp_window **arr, size_t *n,
                     size_t *cap, uint32_t u, uint32_t v, uint32_t w,
                     uint32_t face, double sx, double sy,
                     double t0, double t1, int cross_u, int depth)
{
    double L = edge_len(m, u, v);
    double uw = edge_len(m, u, w), vw = edge_len(m, v, w);
    double x, y2, y, base, len, ang, nsx, nsy;
    double ends[2], vals[2];
    double lo, hi;
    mmp_window nw;
    int i;

    if (!(L > MMP_EPS && uw > MMP_EPS && vw > MMP_EPS)) return 1;
    x = (uw * uw + L * L - vw * vw) / (2.0 * L);
    y2 = uw * uw - x * x;
    y = sqrt(y2 > 0.0 ? y2 : 0.0);
    if (y < MMP_EPS) return 1;   /* degenerate face */

    if (cross_u) {
        base = 0.0;
        len = uw;
        ang = -atan2(y, x);
    } else {
        base = L;
        len = vw;
        ang = -atan2(y, x - L);
    }
    rot(sx - base, sy, ang, &nsx, &nsy);
    /* Canonical frame for `face` on the new edge is the mirror of the
     * rotated frame, so reflect the source back above the edge. */
    nsy = -nsy;
    for (i = 0; i < 2; ++i) {
        double px, py, lam;
        rot((i == 0 ? t0 : t1) - base, 0.0, ang, &px, &py);
        py = -py;   /* same reflection as the source */
        if (fabs(py - nsy) < 1e-14) {
            /* Ray parallel to the target edge: unbounded intersection. */
            vals[i] = (px > nsx) ? 1e18 : -1e18;
            continue;
        }
        lam = -nsy / (py - nsy);
        vals[i] = nsx + lam * (px - nsx);
    }
    lo = vals[0] < vals[1] ? vals[0] : vals[1];
    hi = vals[0] < vals[1] ? vals[1] : vals[0];
    if (lo < 0.0) lo = 0.0;
    if (hi > len) hi = len;
    if (hi - lo < MMP_EPS) return 1;
    ends[0] = lo;
    ends[1] = hi;

    nw.u = cross_u ? u : v;
    nw.v = w;
    nw.face = face;
    nw.t0 = ends[0];
    nw.t1 = ends[1];
    nw.sx = nsx;
    nw.sy = nsy;
    nw.depth = depth + 1;
    return add_window(arr, n, cap, nw);
}

int geomlib_geodesic_mmp(const geomlib_trimesh *mesh,
                         uint32_t source, uint32_t target, double *length)
{
    mmp_window *wins = NULL;
    mmp_edge *edges = NULL;
    size_t nw = 0, cap = 0, ne = 0, i, k, step;
    double best = DBL_MAX;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !length || !mesh->vertices || !mesh->indices ||
        source >= mesh->vertex_count || target >= mesh->vertex_count)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (source == target) { *length = 0.0; return GEOMLIB_MESH_SUCCESS; }

    edges = (mmp_edge *)malloc(mesh->triangle_count * 3 * sizeof(mmp_edge));
    if (!edges) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i = 0; i < mesh->triangle_count; ++i) {
        uint32_t v[3] = {mesh->indices[i * 3], mesh->indices[i * 3 + 1],
                         mesh->indices[i * 3 + 2]};
        for (k = 0; k < 3; ++k) {
            uint32_t a = v[k], b = v[(k + 1) % 3];
            int idx;
            if (a >= mesh->vertex_count || b >= mesh->vertex_count) continue;
            idx = find_edge(edges, ne, a, b);
            if (idx < 0) {
                edges[ne].a = a < b ? a : b;
                edges[ne].b = a < b ? b : a;
                edges[ne].f0 = (uint32_t)i;
                edges[ne].f1 = UINT32_MAX;
                ne++;
            } else if (edges[idx].f1 == UINT32_MAX) {
                edges[idx].f1 = (uint32_t)i;
            }
        }
    }

    /* Seed: unfold every face incident to `source` across its opposite edge. */
    for (i = 0; i < mesh->triangle_count; ++i) {
        uint32_t v[3] = {mesh->indices[i * 3], mesh->indices[i * 3 + 1],
                         mesh->indices[i * 3 + 2]};
        for (k = 0; k < 3; ++k) {
            uint32_t u = v[k], vv = v[(k + 1) % 3], s = v[(k + 2) % 3];
            double L, us, vs, x, y2, y;
            mmp_window wnd;
            if (s != source) continue;
            L = edge_len(mesh, u, vv);
            if (!(L > MMP_EPS)) continue;
            us = edge_len(mesh, u, s);
            vs = edge_len(mesh, vv, s);
            x = (us * us + L * L - vs * vs) / (2.0 * L);
            y2 = us * us - x * x;
            if (y2 <= 0.0) continue;
            y = sqrt(y2);
            wnd.u = u; wnd.v = vv; wnd.face = (uint32_t)i;
            wnd.depth = 0;
            wnd.t0 = 0.0; wnd.t1 = L;
            wnd.sx = x; wnd.sy = y;
            if (!add_window(&wins, &nw, &cap, wnd)) goto cleanup;
        }
    }
    if (!nw) goto cleanup;

    for (step = 0; step < MMP_MAX_WINDOWS && nw <= MMP_MAX_WINDOWS; ++step) {
        size_t bi = (size_t)-1;
        double bk = DBL_MAX;
        mmp_window cur;
        uint32_t opp = UINT32_MAX, w = 0;
        int ei;

        for (i = 0; i < nw; ++i)
            if (!wins[i].processed && wins[i].key < bk) {
                bk = wins[i].key;
                bi = i;
            }
        if (bi == (size_t)-1) break;
        if (bk > best + 1e-9) break;
        wins[bi].processed = 1;
        cur = wins[bi];

        /* Target endpoint check. */
        if (cur.u == target || cur.v == target) {
            double t = (cur.u == target) ? 0.0 :
                       edge_len(mesh, cur.u, cur.v);
            if (t >= cur.t0 - 1e-9 && t <= cur.t1 + 1e-9) {
                double d = dist_at(&cur, t);
                if (d < best) best = d;
            }
        }

        ei = find_edge(edges, ne, cur.u, cur.v);
        if (ei < 0) continue;
        opp = (edges[ei].f0 == cur.face) ? edges[ei].f1 : edges[ei].f0;
        if (opp == UINT32_MAX) continue;   /* boundary edge */
        third_vertex(mesh, opp, cur.u, cur.v, &w);

        /* Crossing the shared edge mirrors the unfolded source. */
        {
            double sx = cur.sx, sy = -cur.sy;
            if (!propagate(mesh, &wins, &nw, &cap, cur.u, cur.v, w, opp,
                           sx, sy, cur.t0, cur.t1, 1, cur.depth)) goto cleanup;
            if (!propagate(mesh, &wins, &nw, &cap, cur.u, cur.v, w, opp,
                           sx, sy, cur.t0, cur.t1, 0, cur.depth)) goto cleanup;
        }
    }

    if (best < DBL_MAX) { *length = best; rc = GEOMLIB_MESH_SUCCESS; }

cleanup:
    free(wins);
    free(edges);
    return rc;
}
