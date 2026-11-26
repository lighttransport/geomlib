/*
 * Convex polygon intersection area: overlap, identical and disjoint.
 */
#include <geomlib/geomlib_polygon_intersection_area.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static void square(geomlib_vec2 out[4], float x, float y, float s)
{
    out[0] = geomlib_vec2_create(x, y);
    out[1] = geomlib_vec2_create(x + s, y);
    out[2] = geomlib_vec2_create(x + s, y + s);
    out[3] = geomlib_vec2_create(x, y + s);
}

int main(void)
{
    geomlib_vec2 a[4], b[4];
    square(a, 0, 0, 1);
    square(b, 0.5f, 0.5f, 1);
    double area = geomlib_polygon_intersection_area(a, 4, b, 4);
    std::printf("intersection area=%.6f expect=0.25\n", area);
    check(std::fabs(area - 0.25) < 1e-5, "half-overlapping squares give 0.25");

    square(b, 0, 0, 1);
    area = geomlib_polygon_intersection_area(a, 4, b, 4);
    check(std::fabs(area - 1.0) < 1e-5, "identical squares give area 1");

    square(b, 3, 3, 1);
    area = geomlib_polygon_intersection_area(a, 4, b, 4);
    check(area == 0.0, "disjoint squares give zero");
    if (failures != 0) {
        std::printf("%d polygon intersection area test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polygon intersection area tests passed\n");
    return 0;
}
