/*
 * geomlib - Triangle / triangle intersection (see geomlib_tri_tri.h)
 *
 * Clean-room separating-plane test with interval overlap on the intersection
 * line, plus a projected 2D test that returns a representative segment in the
 * coplanar case. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_tri_tri.h"

#include <math.h>

#define TT_EPS 1e-7

static void to_arr(geomlib_vec3 v, double o[3])
{
    o[0] = (double)v.x; o[1] = (double)v.y; o[2] = (double)v.z;
}

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double o[3])
{
    o[0] = a[1] * b[2] - a[2] * b[1];
    o[1] = a[2] * b[0] - a[0] * b[2];
    o[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3(const double a[3])
{
    return sqrt(dot3(a, a));
}

static geomlib_tri_tri_result empty_result(void)
{
    geomlib_tri_tri_result r;
    r.intersects = 0;
    r.p0 = r.p1 = geomlib_vec3_create(0, 0, 0);
    return r;
}

static int orient2d(const double a[2], const double b[2], const double c[2])
{
    double d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if (d > TT_EPS) return 1;
    if (d < -TT_EPS) return -1;
    return 0;
}

static int on_segment2d(const double a[2], const double b[2], const double p[2])
{
    return p[0] >= fmin(a[0], b[0]) - TT_EPS &&
           p[0] <= fmax(a[0], b[0]) + TT_EPS &&
           p[1] >= fmin(a[1], b[1]) - TT_EPS &&
           p[1] <= fmax(a[1], b[1]) + TT_EPS;
}

static int point_in_tri2d(const double p[2], double t[3][2])
{
    int d1 = orient2d(t[0], t[1], p);
    int d2 = orient2d(t[1], t[2], p);
    int d3 = orient2d(t[2], t[0], p);
    if (d1 >= 0 && d2 >= 0 && d3 >= 0) return 1;
    if (d1 <= 0 && d2 <= 0 && d3 <= 0) return 1;
    return 0;
}

/* True segment crossing: compute the intersection point when it exists. */
static int seg_point2d(const double a[2], const double b[2],
                       const double c[2], const double d[2], double out[2])
{
    double d1 = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
    double d2 = (b[0]-a[0])*(d[1]-a[1]) - (b[1]-a[1])*(d[0]-a[0]);
    double den;
    if (d1 * d2 >= 0.0 && d1 != 0.0 && d2 != 0.0) return 0;
    den = d1 - d2;
    if (fabs(den) < 1e-18) return 0;
    out[0] = (d1 * d[0] - d2 * c[0]) / den;
    out[1] = (d1 * d[1] - d2 * c[1]) / den;
    if (!on_segment2d(a, b, out) || !on_segment2d(c, d, out)) return 0;
    return 1;
}

static void add_point(double pts[32][2], int *np, const double p[2])
{
    int i;
    for (i = 0; i < *np; ++i)
        if (fabs(pts[i][0] - p[0]) < 1e-5 && fabs(pts[i][1] - p[1]) < 1e-5)
            return;
    if (*np < 32) {
        pts[*np][0] = p[0]; pts[*np][1] = p[1];
        (*np)++;
    }
}

static geomlib_tri_tri_result coplanar_result(double A[3][3],
                                              double B[3][3],
                                              double n[3], double d)
{
    geomlib_tri_tri_result r = empty_result();
    int axis = fabs(n[0]) >= fabs(n[1]) && fabs(n[0]) >= fabs(n[2]) ? 0
             : (fabs(n[1]) >= fabs(n[2]) ? 1 : 2);
    int u = (axis + 1) % 3, v = (axis + 2) % 3;
    double a2[3][2], b2[3][2], pts[32][2];
    int np = 0, i, j, best_i = 0, best_j = 0;
    double best = -1.0;

    for (i = 0; i < 3; ++i) {
        a2[i][0] = A[i][u]; a2[i][1] = A[i][v];
        b2[i][0] = B[i][u]; b2[i][1] = B[i][v];
    }
    for (i = 0; i < 3; ++i) {
        if (point_in_tri2d(a2[i], b2)) add_point(pts, &np, a2[i]);
        if (point_in_tri2d(b2[i], a2)) add_point(pts, &np, b2[i]);
    }
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            double p[2];
            if (seg_point2d(a2[i], a2[(i+1)%3], b2[j], b2[(j+1)%3], p))
                add_point(pts, &np, p);
        }
    if (!np) return r;
    for (i = 0; i < np; ++i)
        for (j = i + 1; j < np; ++j) {
            double dx = pts[i][0] - pts[j][0], dy = pts[i][1] - pts[j][1];
            double dd = dx * dx + dy * dy;
            if (dd > best) { best = dd; best_i = i; best_j = j; }
        }
    {
        double pa[3] = {0, 0, 0}, pb[3] = {0, 0, 0};
        pa[u] = pts[best_i][0]; pa[v] = pts[best_i][1];
        pb[u] = pts[best_j][0]; pb[v] = pts[best_j][1];
        pa[axis] = (d - n[u] * pa[u] - n[v] * pa[v]) / n[axis];
        pb[axis] = (d - n[u] * pb[u] - n[v] * pb[v]) / n[axis];
        r.intersects = 1;
        r.p0 = geomlib_vec3_create((float)pa[0], (float)pa[1], (float)pa[2]);
        r.p1 = geomlib_vec3_create((float)pb[0], (float)pb[1], (float)pb[2]);
    }
    return r;
}

static int interval(double T[3][3], double s[3],
                    double p[3], double dir[3], double dir2,
                    double *tmin, double *tmax)
{
    int i, count = 0;
    for (i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        double den, u, q[3], dv[3], t;
        if (fabs(s[i]) < TT_EPS && fabs(s[j]) < TT_EPS) continue;
        if ((s[i] > TT_EPS && s[j] > TT_EPS) ||
            (s[i] < -TT_EPS && s[j] < -TT_EPS)) continue;
        den = s[i] - s[j];
        if (fabs(den) < 1e-18) continue;
        u = s[i] / den;
        dv[0] = T[j][0] - T[i][0];
        dv[1] = T[j][1] - T[i][1];
        dv[2] = T[j][2] - T[i][2];
        q[0] = T[i][0] + u * dv[0];
        q[1] = T[i][1] + u * dv[1];
        q[2] = T[i][2] + u * dv[2];
        t = ((q[0] - p[0]) * dir[0] + (q[1] - p[1]) * dir[1] +
             (q[2] - p[2]) * dir[2]) / dir2;
        if (count == 0) *tmin = *tmax = t;
        else {
            if (t < *tmin) *tmin = t;
            if (t > *tmax) *tmax = t;
        }
        count++;
    }
    return count;
}

geomlib_tri_tri_result geomlib_triangle_intersection(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2,
    geomlib_vec3 b0, geomlib_vec3 b1, geomlib_vec3 b2)
{
    geomlib_tri_tri_result r = empty_result();
    double A[3][3], B[3][3];
    double e1[3], e2[3], nA[3], nB[3], lA, lB;
    double dA, dB, saB[3], sbA[3];
    double m[3], dir[3], p[3], dir2;
    double tminA, tmaxA, tminB, tmaxB, lo, hi, tol;
    int i, ca, cb;

    to_arr(a0, A[0]); to_arr(a1, A[1]); to_arr(a2, A[2]);
    to_arr(b0, B[0]); to_arr(b1, B[1]); to_arr(b2, B[2]);
    for (i = 0; i < 3; ++i) {
        e1[i] = A[1][i] - A[0][i];
        e2[i] = A[2][i] - A[0][i];
    }
    cross3(e1, e2, nA);
    for (i = 0; i < 3; ++i) {
        e1[i] = B[1][i] - B[0][i];
        e2[i] = B[2][i] - B[0][i];
    }
    cross3(e1, e2, nB);
    lA = norm3(nA); lB = norm3(nB);
    if (!(lA > 1e-18 && lB > 1e-18)) return r;
    for (i = 0; i < 3; ++i) { nA[i] /= lA; nB[i] /= lB; }
    dA = dot3(nA, A[0]);
    dB = dot3(nB, B[0]);
    for (i = 0; i < 3; ++i) {
        saB[i] = dot3(nB, A[i]) - dB;
        sbA[i] = dot3(nA, B[i]) - dA;
    }
    /* Coplanar faces have all signed distances near zero. */
    {
        int coplanar = 1;
        for (i = 0; i < 3; ++i)
            if (fabs(saB[i]) > TT_EPS || fabs(sbA[i]) > TT_EPS) coplanar = 0;
        if (coplanar) return coplanar_result(A, B, nA, dA);
    }
    {
        /* Each triangle must straddle the other's plane: test the two sets
         * independently, otherwise parallel offset planes look straddling. */
        int a_pos = 1, a_neg = 1, b_pos = 1, b_neg = 1;
        for (i = 0; i < 3; ++i) {
            if (saB[i] <= TT_EPS) a_pos = 0;
            if (saB[i] >= -TT_EPS) a_neg = 0;
            if (sbA[i] <= TT_EPS) b_pos = 0;
            if (sbA[i] >= -TT_EPS) b_neg = 0;
        }
        if (a_pos || a_neg || b_pos || b_neg) return r;
    }
    cross3(nA, nB, dir);
    dir2 = dot3(dir, dir);
    if (dir2 < 1e-16) return coplanar_result(A, B, nA, dA);

    for (i = 0; i < 3; ++i) m[i] = dA * nB[i] - dB * nA[i];
    cross3(m, dir, p);
    for (i = 0; i < 3; ++i) p[i] /= dir2;
    ca = interval(A, saB, p, dir, dir2, &tminA, &tmaxA);
    cb = interval(B, sbA, p, dir, dir2, &tminB, &tmaxB);
    if (!ca || !cb) return r;
    lo = fmax(tminA, tminB);
    hi = fmin(tmaxA, tmaxB);
    tol = TT_EPS * (1.0 + fabs(lo) + fabs(hi));
    if (lo > hi + tol) return r;
    r.intersects = 1;
    r.p0 = geomlib_vec3_create((float)(p[0] + lo * dir[0]),
                               (float)(p[1] + lo * dir[1]),
                               (float)(p[2] + lo * dir[2]));
    r.p1 = geomlib_vec3_create((float)(p[0] + hi * dir[0]),
                               (float)(p[1] + hi * dir[1]),
                               (float)(p[2] + hi * dir[2]));
    return r;
}
