/*
 * Circle-line: secant, tangent, miss, and finite segment clipping.
 */
#include <geomlib/geomlib_circle_line.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec2 o = geomlib_vec2_create(0, 0);
    geomlib_circle_line_intersection r = geomlib_circle_line_intersect(
        geomlib_vec2_create(-2, 0), geomlib_vec2_create(2, 0), o, 1.0f);
    std::printf("circle-line kind=%d p0=(%.2f,%.2f) p1=(%.2f,%.2f)\n",
                (int)r.kind, r.p0.x, r.p0.y, r.p1.x, r.p1.y);
    check(r.kind == GEOMLIB_CIRCLE_LINE_TWO_POINTS &&
          std::fabs(r.p0.x + 1.0f) < 1e-5f &&
          std::fabs(r.p1.x - 1.0f) < 1e-5f, "secant gives ±1");

    r = geomlib_circle_line_intersect(
        geomlib_vec2_create(-2, 1), geomlib_vec2_create(2, 1), o, 1.0f);
    check(r.kind == GEOMLIB_CIRCLE_LINE_TANGENT &&
          std::fabs(r.p0.x) < 1e-5f && std::fabs(r.p0.y - 1.0f) < 1e-5f,
          "tangent line touches at (0,1)");

    r = geomlib_circle_line_intersect(
        geomlib_vec2_create(-2, 2), geomlib_vec2_create(2, 2), o, 1.0f);
    check(r.kind == GEOMLIB_CIRCLE_LINE_NONE, "distant line misses");

    /* Segment [0.5,2] along x contains only x=1. */
    r = geomlib_circle_segment_intersect(
        geomlib_vec2_create(0.5f, 0), geomlib_vec2_create(2, 0), o, 1.0f);
    std::printf("circle-segment kind=%d p0.x=%.2f\n", (int)r.kind, r.p0.x);
    check(r.kind == GEOMLIB_CIRCLE_LINE_TANGENT &&
          std::fabs(r.p0.x - 1.0f) < 1e-5f,
          "segment clipping keeps only in-range intersections");
    if (failures != 0) {
        std::printf("%d circle-line test(s) failed\n", failures);
        return 1;
    }
    std::printf("All circle-line tests passed\n");
    return 0;
}
