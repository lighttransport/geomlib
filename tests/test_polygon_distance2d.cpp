/*
 * 2D polygon distance: separated by a gap, overlapping, and diagonal.
 */
#include <geomlib/geomlib_polygon_distance2d.h>

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
    geomlib_polygon_distance2d_result r;
    square(a, 0, 0, 1);
    square(b, 2, 0, 1);

    int ok = geomlib_polygon_distance2d(a, 4, b, 4, &r);
    std::printf("poly dist=%.4f points=(%.2f,%.2f)-(%.2f,%.2f)\n",
                r.distance, r.point_a.x, r.point_a.y, r.point_b.x, r.point_b.y);
    check(ok && std::fabs(r.distance - 1.0f) < 1e-4f &&
          std::fabs(r.point_a.x - 1.0f) < 1e-4f &&
          std::fabs(r.point_b.x - 2.0f) < 1e-4f,
          "unit gap between squares with witnesses");

    square(b, 0.5f, 0.5f, 1);
    ok = geomlib_polygon_distance2d(a, 4, b, 4, &r);
    check(ok && r.distance == 0.0f, "overlapping polygons have distance zero");

    square(b, 2, 2, 1);
    ok = geomlib_polygon_distance2d(a, 4, b, 4, &r);
    std::printf("diagonal poly dist=%.6f\n", r.distance);
    check(ok && std::fabs(r.distance - std::sqrt(2.0)) < 1e-4f,
          "diagonal gap is sqrt(2)");
    if (failures != 0) {
        std::printf("%d polygon distance test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polygon distance tests passed\n");
    return 0;
}
