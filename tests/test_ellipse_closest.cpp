/*
 * Closest point on ellipse: axial rays, interior query and a rotated ellipse.
 */
#include <geomlib/geomlib_ellipse_closest.h>

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
    geomlib_ellipse_closest r;

    int ok = geomlib_ellipse_closest_point(
        geomlib_vec2_create(3, 0), o, 2.0f, 1.0f, 0.0f, &r);
    std::printf("closest to (3,0)=(%.3f,%.3f) d=%.4f\n",
                r.point.x, r.point.y, r.distance);
    check(ok && std::fabs(r.point.x - 2.0f) < 1e-3f &&
          std::fabs(r.distance - 1.0f) < 1e-3f,
          "point beyond the major axis projects to (2,0)");

    ok = geomlib_ellipse_closest_point(
        geomlib_vec2_create(0, 2), o, 2.0f, 1.0f, 0.0f, &r);
    check(ok && std::fabs(r.point.y - 1.0f) < 1e-3f &&
          std::fabs(r.distance - 1.0f) < 1e-3f,
          "point beyond the minor axis projects to (0,1)");

    ok = geomlib_ellipse_closest_point(o, o, 2.0f, 1.0f, 0.0f, &r);
    check(ok && std::fabs(r.distance - 1.0f) < 1e-3f,
          "interior centre projects to the nearer minor-axis point");

    float ang = (float)(3.14159265358979323846 / 4.0);
    float s = (float)(std::sqrt(2.0) / 2.0);
    ok = geomlib_ellipse_closest_point(
        geomlib_vec2_create(3 * s, 3 * s), o, 2.0f, 1.0f, ang, &r);
    std::printf("rotated closest=(%.3f,%.3f) d=%.4f\n",
                r.point.x, r.point.y, r.distance);
    check(ok && std::fabs(r.distance - 1.0f) < 1e-3f &&
          std::fabs(r.point.x - 2 * s) < 1e-3f,
          "rotated major-axis projection is exact");
    if (failures != 0) {
        std::printf("%d ellipse closest test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ellipse closest tests passed\n");
    return 0;
}
