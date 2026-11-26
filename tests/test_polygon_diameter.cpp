/*
 * Convex diameter: square, rectangle and a dense circle approximation.
 */
#include <geomlib/geomlib_polygon_diameter.h>

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
    double d = geomlib_convex_polygon_diameter(sq, 4);
    std::printf("square diameter=%.6f expect=%.6f\n", d, std::sqrt(2.0));
    check(std::fabs(d - std::sqrt(2.0)) < 1e-5, "unit square diagonal");

    geomlib_vec2 rect[4] = {
        geomlib_vec2_create(0, 0), geomlib_vec2_create(4, 0),
        geomlib_vec2_create(4, 2), geomlib_vec2_create(0, 2)
    };
    d = geomlib_convex_polygon_diameter(rect, 4);
    check(std::fabs(d - std::sqrt(20.0)) < 1e-5, "4x2 rectangle diagonal");

    std::vector<geomlib_vec2> circle;
    for (int i = 0; i < 128; ++i) {
        double a = 2.0 * 3.14159265358979323846 * i / 128.0;
        circle.push_back(geomlib_vec2_create(
            (float)(2.0 * cos(a)), (float)(2.0 * sin(a))));
    }
    d = geomlib_convex_polygon_diameter(circle.data(), circle.size());
    std::printf("circle diameter=%.6f expect=4\n", d);
    check(std::fabs(d - 4.0) < 0.02, "circle diameter is about 4");
    if (failures != 0) {
        std::printf("%d polygon diameter test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polygon diameter tests passed\n");
    return 0;
}
