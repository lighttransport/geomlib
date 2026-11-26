/*
 * Polyline closest point: interior projection, endpoint clamping and a
 * single-point polyline.
 */
#include <geomlib/geomlib_polyline_closest.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 pts[3] = {
        geomlib_vec3_create(0, 0, 0),
        geomlib_vec3_create(2, 0, 0),
        geomlib_vec3_create(2, 2, 0)
    };
    geomlib_polyline_closest r;

    int ok = geomlib_polyline_closest_point(
        pts, 3, geomlib_vec3_create(1, 1, 0), &r);
    std::printf("closest point=(%.2f,%.2f,%.2f) d=%.4f seg=%zu t=%.2f\n",
                r.point.x, r.point.y, r.point.z, r.distance, r.segment, r.t);
    check(ok && std::fabs(r.point.x - 1.0f) < 1e-5f &&
          std::fabs(r.point.y) < 1e-5f && std::fabs(r.distance - 1.0f) < 1e-5f &&
          r.segment == 0 && std::fabs(r.t - 0.5f) < 1e-5f,
          "interior projection onto the first segment");

    ok = geomlib_polyline_closest_point(
        pts, 3, geomlib_vec3_create(3, 2, 0), &r);
    std::printf("clamped point=(%.2f,%.2f) seg=%zu t=%.2f\n",
                r.point.x, r.point.y, r.segment, r.t);
    check(ok && std::fabs(r.point.x - 2.0f) < 1e-5f &&
          std::fabs(r.point.y - 2.0f) < 1e-5f &&
          std::fabs(r.distance - 1.0f) < 1e-5f &&
          r.segment == 1 && std::fabs(r.t - 1.0f) < 1e-5f,
          "endpoint clamping on the second segment");

    ok = geomlib_polyline_closest_point(
        pts, 1, geomlib_vec3_create(3, 4, 0), &r);
    check(ok && std::fabs(r.distance - 5.0f) < 1e-5f,
          "single-point polyline distance");
    if (failures != 0) {
        std::printf("%d polyline closest test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polyline closest tests passed\n");
    return 0;
}
