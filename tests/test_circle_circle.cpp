/*
 * Circle-circle: analytic two-point intersection, tangency, disjointness and
 * containment.
 */
#include <geomlib/geomlib_circle_circle.h>

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
    geomlib_vec2 one = geomlib_vec2_create(1, 0);

    geomlib_circle2d_intersection r = geomlib_circle2d_intersect(o, 1.0f, one, 1.0f);
    std::printf("circle-circle kind=%d p0=(%.3f,%.3f)\n",
                (int)r.kind, r.p0.x, r.p0.y);
    check(r.kind == GEOMLIB_CIRCLE2D_TWO_POINTS &&
          std::fabs(r.p0.x - 0.5f) < 1e-4f &&
          std::fabs(std::fabs(r.p0.y) - std::sqrt(3.0) / 2.0) < 1e-4f,
          "unit circles one apart intersect analytically");

    r = geomlib_circle2d_intersect(o, 1.0f, geomlib_vec2_create(2, 0), 1.0f);
    check(r.kind == GEOMLIB_CIRCLE2D_TANGENT &&
          std::fabs(r.p0.x - 1.0f) < 1e-4f &&
          std::fabs(r.p0.y) < 1e-4f, "externally tangent at (1,0)");

    r = geomlib_circle2d_intersect(o, 1.0f, geomlib_vec2_create(3, 0), 1.0f);
    check(r.kind == GEOMLIB_CIRCLE2D_NONE, "disjoint circles");

    r = geomlib_circle2d_intersect(o, 3.0f, geomlib_vec2_create(0.5f, 0), 1.0f);
    check(r.kind == GEOMLIB_CIRCLE2D_CONTAINED, "one circle inside another");

    r = geomlib_circle2d_intersect(o, 2.0f, o, 2.0f);
    check(r.kind == GEOMLIB_CIRCLE2D_COINCIDENT, "coincident circles");
    if (failures != 0) {
        std::printf("%d circle-circle test(s) failed\n", failures);
        return 1;
    }
    std::printf("All circle-circle tests passed\n");
    return 0;
}
