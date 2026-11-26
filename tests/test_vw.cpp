/*
 * Visvalingam-Whyatt: a large spike survives while collinear points are
 * removed; effective areas are exact.
 */
#include <geomlib/geomlib_vw.h>

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
    std::vector<geomlib_vec3> pts;
    for (int i = 0; i < 8; ++i)
        pts.push_back(geomlib_vec3_create((float)i, (i == 4) ? 5.0f : 0.0f, 0));
    double a = geomlib_vw_effective_area(pts.data(), pts.size(), 4);
    std::printf("vw effective area at spike=%.4f expect=5\n", a);
    check(std::fabs(a - 5.0) < 1e-6, "effective area is exact");

    std::vector<uint32_t> idx(pts.size());
    size_t kept = 0;
    int rc = geomlib_simplify_polyline_vw(pts.data(), pts.size(), 0.1,
                                          idx.data(), &kept);
    std::printf("vw rc=%d kept=%zu indices=%u,%u,%u\n", rc, kept,
                idx[0], idx[1], idx[2]);
    check(rc == GEOMLIB_MESH_SUCCESS && kept == 5 &&
          idx[0] == 0 && idx[1] == 3 && idx[2] == 4 &&
          idx[3] == 5 && idx[4] == 7,
          "collinear points removed, spike and its ramps survive");

    std::vector<geomlib_vec3> line;
    for (int i = 0; i < 10; ++i)
        line.push_back(geomlib_vec3_create((float)i, 0, 0));
    rc = geomlib_simplify_polyline_vw(line.data(), line.size(), 0.0,
                                      idx.data(), &kept);
    check(rc == GEOMLIB_MESH_SUCCESS && kept == 2 && idx[0] == 0 && idx[1] == 9,
          "straight line keeps only endpoints");
    if (failures != 0) {
        std::printf("%d VW test(s) failed\n", failures);
        return 1;
    }
    std::printf("All VW tests passed\n");
    return 0;
}
