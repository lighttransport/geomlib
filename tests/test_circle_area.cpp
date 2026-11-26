/*
 * Circle lens area: containment, tangency, partial overlap and union.
 */
#include <geomlib/geomlib_circle_area.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    double pi = std::acos(-1.0);
    double a = geomlib_circle_intersection_area(1, 1, 0);
    std::printf("coincident lens=%.6f expect=%.6f\n", a, pi);
    check(std::fabs(a - pi) < 1e-9, "coincident unit circles give pi");

    check(geomlib_circle_intersection_area(1, 1, 2) == 0.0,
          "externally tangent circles have zero overlap");
    check(geomlib_circle_intersection_area(1, 1, 3) == 0.0,
          "disjoint circles have zero overlap");

    a = geomlib_circle_intersection_area(1, 1, 1);
    std::printf("unit lens at d=1=%.6f\n", a);
    check(std::fabs(a - (2.0 * pi / 3.0 - std::sqrt(3.0) / 2.0)) < 1e-9,
          "lens area at d=1 is 2pi/3 - sqrt(3)/2");

    a = geomlib_circle_union_area(1, 1, 1);
    check(std::fabs(a - (2.0 * pi - (2.0 * pi / 3.0 - std::sqrt(3.0) / 2.0))) < 1e-9,
          "union area is the sum minus the lens");
    if (failures != 0) {
        std::printf("%d circle area test(s) failed\n", failures);
        return 1;
    }
    std::printf("All circle area tests passed\n");
    return 0;
}
