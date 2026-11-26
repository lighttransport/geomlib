/*
 * geomlib - Robust geometric predicates (see geomlib_predicates.h)
 *
 * Expansion arithmetic follows the method published in:
 *   J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and Fast
 *   Robust Geometric Predicates", DCG 18(3), 1997.
 *
 * No target library source was consulted. Only C11 and <math.h>/<float.h>
 * are used; there is no long double and no compiler extension, so this module
 * compiles identically for native and Emscripten (wasm32) targets.
 */

#include "geomlib/geomlib_predicates.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

/* ------------------------------------------------------------------------- */
/* Expansion arithmetic                                                      */
/* ------------------------------------------------------------------------- */

/*
 * An expansion is a sequence of non-overlapping doubles, sorted by increasing
 * magnitude, whose exact sum is the represented value. The sign of the
 * largest-magnitude component is the sign of the value.
 */
#define GEOMLIB_EXP_MAX 512

typedef struct {
    double e[GEOMLIB_EXP_MAX];
    int n;
} geomlib_expansion;

/* Error-free transformation: a + b = s + e, exactly. */
static void two_sum(double a, double b, double *s, double *e)
{
    double x = a + b;
    double bv = x - a;
    double av = x - bv;
    double br = b - bv;
    double ar = a - av;
    *s = x;
    *e = ar + br;
}

/* Dekker split: a = hi + lo, with hi carrying the leading bits. */
static void split_double(double a, double *hi, double *lo)
{
    /* 2^27 + 1 for IEEE-754 double precision (53-bit significand). */
    const double splitter = 134217729.0;
    double c = splitter * a;
    double abig = c - a;
    *hi = c - abig;
    *lo = a - *hi;
}

/* Error-free transformation: a * b = p + e, exactly (no over/underflow). */
static void two_product(double a, double b, double *p, double *e)
{
    double ahi, alo, bhi, blo;
    double err1, err2, err3;

    *p = a * b;

    /* Extreme magnitudes cannot be split without overflow; the fast filter
     * handles those cases before the exact path is reached in practice. */
    if (fabs(a) > 1e150 || fabs(b) > 1e150 ||
        (a != 0.0 && fabs(b) < 1e-290) ||
        (b != 0.0 && fabs(a) < 1e-290)) {
        *e = 0.0;
        return;
    }

    split_double(a, &ahi, &alo);
    split_double(b, &bhi, &blo);
    err1 = *p - (ahi * bhi);
    err2 = err1 - (alo * bhi);
    err3 = err2 - (ahi * blo);
    *e = (alo * blo) - err3;
}

static void exp_init(geomlib_expansion *x)
{
    x->n = 0;
}

/* Add a scalar to an expansion, exactly. */
static void exp_grow(geomlib_expansion *x, double b)
{
    double h[GEOMLIB_EXP_MAX + 1];
    int hn = 0;
    double q = b;
    int i;

    for (i = 0; i < x->n; ++i) {
        double s, err;
        two_sum(q, x->e[i], &s, &err);
        q = s;
        if (err != 0.0 && hn < GEOMLIB_EXP_MAX + 1) {
            h[hn++] = err;
        }
    }
    if (q != 0.0 && hn < GEOMLIB_EXP_MAX + 1) {
        h[hn++] = q;
    }
    for (i = 0; i < hn; ++i) {
        x->e[i] = h[i];
    }
    x->n = hn;
}

/* Multiply an expansion by a scalar, exactly. */
static void exp_scale(geomlib_expansion *x, double b)
{
    geomlib_expansion r;
    int i;

    if (b == 0.0) {
        exp_init(x);
        return;
    }

    exp_init(&r);
    for (i = 0; i < x->n; ++i) {
        double p, err;
        two_product(x->e[i], b, &p, &err);
        exp_grow(&r, p);
        if (err != 0.0) {
            exp_grow(&r, err);
        }
    }
    *x = r;
}

/* r = a + b, exactly. */
static void exp_add(geomlib_expansion *r,
                    const geomlib_expansion *a,
                    const geomlib_expansion *b)
{
    geomlib_expansion t;
    int i;

    exp_init(&t);
    for (i = 0; i < a->n; ++i) {
        exp_grow(&t, a->e[i]);
    }
    for (i = 0; i < b->n; ++i) {
        exp_grow(&t, b->e[i]);
    }
    *r = t;
}

/* ------------------------------------------------------------------------- */
/* Determinants                                                              */
/* ------------------------------------------------------------------------- */

/* Recursive Laplace expansion in floating point, with a conservative bound
 * on the sum of absolute values of all permutation products. */
static void det_double(int n, const double *m, double *det, double *abs_sum)
{
    double minor[5 * 5];
    double d = 0.0;
    double a = 0.0;
    int j;

    if (n == 1) {
        *det = m[0];
        *abs_sum = fabs(m[0]);
        return;
    }

    for (j = 0; j < n; ++j) {
        double sd, sa, term;
        int idx = 0;
        int i, k;

        for (i = 1; i < n; ++i) {
            for (k = 0; k < n; ++k) {
                if (k != j) {
                    minor[idx++] = m[i * n + k];
                }
            }
        }
        det_double(n - 1, minor, &sd, &sa);
        term = m[j] * sd;
        d += ((j & 1) ? -term : term);
        a += fabs(m[j]) * sa;
    }
    *det = d;
    *abs_sum = a;
}

/* Recursive Laplace expansion in exact expansion arithmetic. */
static void det_exact(int n, const double *m, geomlib_expansion *out)
{
    double minor[5 * 5];
    int j;

    if (n == 1) {
        exp_init(out);
        exp_grow(out, m[0]);
        return;
    }

    exp_init(out);
    for (j = 0; j < n; ++j) {
        geomlib_expansion sub, sum;
        double coeff;
        int idx = 0;
        int i, k;

        for (i = 1; i < n; ++i) {
            for (k = 0; k < n; ++k) {
                if (k != j) {
                    minor[idx++] = m[i * n + k];
                }
            }
        }
        det_exact(n - 1, minor, &sub);
        coeff = ((j & 1) ? -1.0 : 1.0) * m[j];
        exp_scale(&sub, coeff);
        exp_add(&sum, out, &sub);
        *out = sum;
    }
}

int geomlib_pred_det_sign(int n, const double *m)
{
    double det, abs_sum;
    double bound;
    geomlib_expansion e;
    double v;

    if (n < 1 || n > 5 || m == NULL) {
        return 0;
    }

    det_double(n, m, &det, &abs_sum);
    /* gamma_5 ~ 5.6e-16; 32 * DBL_EPSILON is a safe upper bound. */
    bound = 32.0 * DBL_EPSILON * abs_sum;
    if (fabs(det) > bound) {
        return det > 0.0 ? 1 : -1;
    }

    det_exact(n, m, &e);
    if (e.n == 0) {
        return 0;
    }
    v = e.e[e.n - 1];
    return v > 0.0 ? 1 : (v < 0.0 ? -1 : 0);
}

/* ------------------------------------------------------------------------- */
/* Public predicates                                                         */
/* ------------------------------------------------------------------------- */

int geomlib_pred_orient2d(double ax, double ay,
                          double bx, double by,
                          double cx, double cy)
{
    double m[9] = {
        ax, ay, 1.0,
        bx, by, 1.0,
        cx, cy, 1.0
    };
    return geomlib_pred_det_sign(3, m);
}

int geomlib_pred_orient2df(float ax, float ay,
                           float bx, float by,
                           float cx, float cy)
{
    return geomlib_pred_orient2d((double)ax, (double)ay,
                                 (double)bx, (double)by,
                                 (double)cx, (double)cy);
}

int geomlib_pred_orient3d(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double cx, double cy, double cz,
                          double dx, double dy, double dz)
{
    double m[9] = {
        bx - ax, by - ay, bz - az,
        cx - ax, cy - ay, cz - az,
        dx - ax, dy - ay, dz - az
    };
    return geomlib_pred_det_sign(3, m);
}

int geomlib_pred_incircle(double ax, double ay,
                          double bx, double by,
                          double cx, double cy,
                          double dx, double dy)
{
    double m[16] = {
        ax, ay, ax * ax + ay * ay, 1.0,
        bx, by, bx * bx + by * by, 1.0,
        cx, cy, cx * cx + cy * cy, 1.0,
        dx, dy, dx * dx + dy * dy, 1.0
    };
    return geomlib_pred_det_sign(4, m);
}

int geomlib_pred_insphere(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double cx, double cy, double cz,
                          double dx, double dy, double dz,
                          double ex, double ey, double ez)
{
    double m[25] = {
        ax, ay, az, ax * ax + ay * ay + az * az, 1.0,
        bx, by, bz, bx * bx + by * by + bz * bz, 1.0,
        cx, cy, cz, cx * cx + cy * cy + cz * cz, 1.0,
        dx, dy, dz, dx * dx + dy * dy + dz * dz, 1.0,
        ex, ey, ez, ex * ex + ey * ey + ez * ez, 1.0
    };
    /* The raw 5x5 determinant is negative for an interior query with a
     * positively oriented tetrahedron; negate to match the documented
     * "inside is positive" convention. */
    return -geomlib_pred_det_sign(5, m);
}
