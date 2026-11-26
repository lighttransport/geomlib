/*
 * Robust predicate tests.
 *
 * Golden signs are computed independently with Python's fractions.Fraction
 * (see scripts/gen_pred_vectors.py -> tests/generated/pred_vectors.h), so the
 * adaptive exact path is validated against exact rational arithmetic rather
 * than against another floating-point implementation.
 */
#include <geomlib/geomlib_predicates.h>

#include <cstdio>

#include "generated/pred_vectors.h"

static int failures = 0;

static void check(int got, int expected, const char *name, int index)
{
    if (got != expected) {
        std::printf("FAIL %s[%d]: got %d expected %d\n", name, index, got, expected);
        failures++;
    }
}

int main(void)
{
    /* Basic orientation conventions. */
    check(geomlib_pred_orient2d(0, 0, 1, 0, 0, 1), 1, "orient2d-ccw", 0);
    check(geomlib_pred_orient2d(0, 0, 0, 1, 1, 0), -1, "orient2d-cw", 0);
    check(geomlib_pred_orient2d(0, 0, 1, 1, 2, 2), 0, "orient2d-line", 0);
    check(geomlib_pred_orient2df(0.f, 0.f, 1.f, 0.f, 0.f, 1.f), 1, "orient2df", 0);

    check(geomlib_pred_orient3d(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1), 1, "orient3d-tet", 0);
    check(geomlib_pred_orient3d(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0), 0, "orient3d-flat", 0);

    check(geomlib_pred_incircle(0, 0, 1, 0, 0, 1, 0.2, 0.2), 1, "incircle-in", 0);
    check(geomlib_pred_incircle(0, 0, 1, 0, 0, 1, 5, 5), -1, "incircle-out", 0);
    check(geomlib_pred_incircle(0, 0, 1, 0, 0, 1, 1, 1), 0, "incircle-cocirc", 0);

    check(geomlib_pred_insphere(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                                0.1, 0.1, 0.1), 1, "insphere-in", 0);
    check(geomlib_pred_insphere(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                                5, 5, 5), -1, "insphere-out", 0);

    /* Exact-rational golden vectors (near-degenerate integer inputs). */
    for (int i = 0; i < COUNT_o2; ++i) {
        const double *v = vec_o2[i];
        check(geomlib_pred_orient2d(v[0], v[1], v[2], v[3], v[4], v[5]),
              expect_o2[i], "o2", i);
    }
    for (int i = 0; i < COUNT_o3; ++i) {
        const double *v = vec_o3[i];
        check(geomlib_pred_orient3d(v[0],v[1],v[2], v[3],v[4],v[5],
                                    v[6],v[7],v[8], v[9],v[10],v[11]),
              expect_o3[i], "o3", i);
    }
    for (int i = 0; i < COUNT_ic; ++i) {
        const double *v = vec_ic[i];
        check(geomlib_pred_incircle(v[0],v[1], v[2],v[3], v[4],v[5], v[6],v[7]),
              expect_ic[i], "ic", i);
    }
    for (int i = 0; i < COUNT_is; ++i) {
        const double *v = vec_is[i];
        check(geomlib_pred_insphere(v[0],v[1],v[2], v[3],v[4],v[5],
                                    v[6],v[7],v[8], v[9],v[10],v[11],
                                    v[12],v[13],v[14]),
              expect_is[i], "is", i);
    }

    /* Generic determinant sign agrees with the dedicated predicates. */
    double m[9] = {0,0,1, 1,0,1, 0,1,1};
    check(geomlib_pred_det_sign(3, m), 1, "det_sign", 0);
    check(geomlib_pred_det_sign(0, m), 0, "det_sign-invalid", 0);

    if (failures != 0) {
        std::printf("%d predicate test(s) failed\n", failures);
        return 1;
    }
    std::printf("All robust predicate tests passed\n");
    return 0;
}
