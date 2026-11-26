/*
 * Ellipse: axis-aligned implicit values and containment under rotation.
 */
#include <geomlib/geomlib_ellipse.h>

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
    double v = geomlib_ellipse_implicit(
        geomlib_vec2_create(1, 0), o, 2.0f, 1.0f, 0.0f);
    std::printf("ellipse implicit at (1,0)=%.4f\n", v);
    check(std::fabs(v - 0.25) < 1e-9, "axis-aligned implicit value");
    check(geomlib_point_in_ellipse(geomlib_vec2_create(1, 0), o, 2, 1, 0),
          "interior point");
    check(!geomlib_point_in_ellipse(geomlib_vec2_create(0, 2), o, 2, 1, 0),
          "point beyond the minor axis is outside");
    check(geomlib_point_in_ellipse(geomlib_vec2_create(2, 0), o, 2, 1, 0),
          "boundary point is inside");

    float a = (float)(3.14159265358979323846 / 4.0);
    float s = (float)(std::sqrt(2.0) / 2.0);
    check(geomlib_point_in_ellipse(
              geomlib_vec2_create(1.5f * s, 1.5f * s), o, 2, 1, a),
          "point on the rotated major axis is inside");
    check(!geomlib_point_in_ellipse(
              geomlib_vec2_create(1.5f * s, -1.5f * s), o, 2, 1, a),
          "point on the rotated minor axis is outside");
    if (failures != 0) {
        std::printf("%d ellipse test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ellipse tests passed\n");
    return 0;
}
