/*
 * Marching squares test: a circle implicit field yields a closed contour
 * whose length approximates 2*pi*r.
 */
#include <geomlib/geomlib_marching_squares.h>

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
    const size_t n = 64;
    const float h = 4.0f/n;
    std::vector<float> field(n*n);
    for (size_t y = 0; y < n; ++y)
        for (size_t x = 0; x < n; ++x) {
            float px = -2.0f + (float)x*h;
            float py = -2.0f + (float)y*h;
            field[x+n*y] = px*px + py*py - 1.0f;
        }
    geomlib_contour_result r = geomlib_contour_result_create();
    int rc = geomlib_marching_squares(field.data(), n, n, h,
                                      geomlib_vec2_create(-2.0f,-2.0f), &r);
    double length = 0.0;
    for (size_t p = 0; p < r.polyline_count; ++p)
        for (uint32_t i = 1; i < r.counts[p]; ++i) {
            geomlib_vec2 a = r.points[r.starts[p]+i-1];
            geomlib_vec2 b = r.points[r.starts[p]+i];
            length += std::sqrt((b.x-a.x)*(b.x-a.x) + (b.y-a.y)*(b.y-a.y));
        }
    std::printf("marching squares rc=%d contours=%zu points=%zu length=%.4f expected=%.4f\n",
                rc, r.polyline_count, r.point_count, length, 6.2831853);
    check(rc == GEOMLIB_MESH_SUCCESS, "marching squares succeeds");
    check(r.polyline_count >= 1 && r.point_count > 20, "contour has points");
    check(std::fabs(length - 6.2831853) < 0.2, "circle perimeter approximates 2*pi");
    geomlib_contour_result_free(&r);
    if (failures != 0) {
        std::printf("%d marching squares test(s) failed\n", failures);
        return 1;
    }
    std::printf("All marching squares tests passed\n");
    return 0;
}
