/*
 * 2D line intersection: crossing, parallel and collinear.
 */
#include <geomlib/geomlib_line2d.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_line2d_intersection r = geomlib_line2d_intersect(
        geomlib_vec2_create(0, 0), geomlib_vec2_create(2, 0),
        geomlib_vec2_create(1, -1), geomlib_vec2_create(1, 1));
    std::printf("line2d kind=%d point=(%.2f,%.2f) t=%.2f u=%.2f\n",
                (int)r.kind, r.point.x, r.point.y, r.t, r.u);
    check(r.kind == GEOMLIB_LINE2D_POINT &&
          std::fabs(r.point.x - 1.0f) < 1e-5f &&
          std::fabs(r.point.y) < 1e-5f &&
          std::fabs(r.t - 0.5) < 1e-9 && std::fabs(r.u - 0.5) < 1e-9,
          "crossing lines give exact point and parameters");

    r = geomlib_line2d_intersect(
        geomlib_vec2_create(0, 0), geomlib_vec2_create(2, 0),
        geomlib_vec2_create(0, 1), geomlib_vec2_create(2, 1));
    check(r.kind == GEOMLIB_LINE2D_NONE, "parallel distinct lines miss");

    r = geomlib_line2d_intersect(
        geomlib_vec2_create(0, 0), geomlib_vec2_create(2, 0),
        geomlib_vec2_create(1, 0), geomlib_vec2_create(3, 0));
    check(r.kind == GEOMLIB_LINE2D_COLLINEAR, "collinear lines are classified");
    if (failures != 0) {
        std::printf("%d line2d test(s) failed\n", failures);
        return 1;
    }
    std::printf("All line2d tests passed\n");
    return 0;
}
