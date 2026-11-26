/*
 * Ellipse measures: circle limits and Ramanujan perimeter accuracy.
 */
#include <geomlib/geomlib_ellipse_measures.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    double p = geomlib_ellipse_perimeter(1.0f, 1.0f);
    std::printf("circle perimeter=%.8f expect=%.8f\n",
                p, 2.0 * std::acos(-1.0));
    check(std::fabs(p - 2.0 * std::acos(-1.0)) < 1e-9,
          "circle perimeter is exact");

    p = geomlib_ellipse_perimeter(2.0f, 1.0f);
    std::printf("ellipse perimeter=%.6f expect~9.68845\n", p);
    check(std::fabs(p - 9.68845) < 1e-4, "Ramanujan approximation accuracy");

    double a = geomlib_ellipse_area(2.0f, 1.0f);
    check(std::fabs(a - 2.0 * std::acos(-1.0)) < 1e-9, "area is pi*a*b");
    if (failures != 0) {
        std::printf("%d ellipse measure test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ellipse measure tests passed\n");
    return 0;
}
