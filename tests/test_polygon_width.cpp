/*
 * Convex width: unit square, rectangle, equilateral triangle and circle.
 */
#include <geomlib/geomlib_polygon_width.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec2 sq[4] = {
        geomlib_vec2_create(0, 0), geomlib_vec2_create(1, 0),
        geomlib_vec2_create(1, 1), geomlib_vec2_create(0, 1)
    };
    double w = geomlib_convex_polygon_min_width(sq, 4);
    std::printf("square width=%.6f expect=1\n", w);
    check(std::fabs(w - 1.0) < 1e-5, "unit square width is 1");

    geomlib_vec2 rect[4] = {
        geomlib_vec2_create(0, 0), geomlib_vec2_create(4, 0),
        geomlib_vec2_create(4, 2), geomlib_vec2_create(0, 2)
    };
    w = geomlib_convex_polygon_min_width(rect, 4);
    check(std::fabs(w - 2.0) < 1e-5, "4x2 rectangle width is 2");

    geomlib_vec2 tri[3] = {
        geomlib_vec2_create(0, 0), geomlib_vec2_create(2, 0),
        geomlib_vec2_create(1, (float)std::sqrt(3.0))
    };
    w = geomlib_convex_polygon_min_width(tri, 3);
    std::printf("triangle width=%.6f expect=%.6f\n", w, std::sqrt(3.0));
    check(std::fabs(w - std::sqrt(3.0)) < 1e-4, "equilateral triangle height");

    std::vector<geomlib_vec2> circle;
    for (int i = 0; i < 128; ++i) {
        double a = 2.0 * 3.14159265358979323846 * i / 128.0;
        circle.push_back(geomlib_vec2_create(
            (float)(2.0 * cos(a)), (float)(2.0 * sin(a))));
    }
    w = geomlib_convex_polygon_min_width(circle.data(), circle.size());
    std::printf("circle width=%.6f expect about 4\n", w);
    check(std::fabs(w - 4.0) < 0.02, "circle width is about 4");
    if (failures != 0) {
        std::printf("%d polygon width test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polygon width tests passed\n");
    return 0;
}
