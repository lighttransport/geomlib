/*
 * Douglas-Peucker tests: collinear points collapse to endpoints, zigzags keep
 * peaks above tolerance, endpoints are preserved and error is bounded.
 */
#include <geomlib/geomlib_simplify.h>

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
    std::vector<geomlib_vec3> line;
    for (int i = 0; i < 50; ++i)
        line.push_back(geomlib_vec3_create((float)i, 0.0f, 0.0f));
    std::vector<uint32_t> idx(line.size());
    size_t n = 0;
    int rc = geomlib_simplify_polyline(line.data(), line.size(), 0.01f,
                                       idx.data(), &n);
    std::printf("simplify line rc=%d kept=%zu\n", rc, n);
    check(rc == GEOMLIB_MESH_SUCCESS && n == 2, "collinear line keeps two endpoints");
    check(idx[0] == 0 && idx[1] == line.size()-1, "endpoints preserved");
    check(geomlib_simplify_max_error(line.data(), line.size(), idx.data(), n) < 1e-5f,
          "line error is zero");

    std::vector<geomlib_vec3> zig;
    for (int i = 0; i < 21; ++i)
        zig.push_back(geomlib_vec3_create((float)i,
                                          (i%5==2) ? 0.5f : 0.01f, 0.0f));
    rc = geomlib_simplify_polyline(zig.data(), zig.size(), 0.1f, idx.data(), &n);
    float err = geomlib_simplify_max_error(zig.data(), zig.size(), idx.data(), n);
    std::printf("simplify zigzag rc=%d kept=%zu error=%.4f\n", rc, n, err);
    check(rc == GEOMLIB_MESH_SUCCESS, "zigzag simplifies");
    check(n > 2 && n < zig.size(), "interior peaks retained, others removed");
    check(idx[0] == 0 && idx[n-1] == zig.size()-1, "zigzag endpoints preserved");
    check(err <= 0.1001f, "max deviation within tolerance");

    std::vector<geomlib_vec3> tri = {
        geomlib_vec3_create(0,0,0), geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0)
    };
    rc = geomlib_simplify_polyline(tri.data(), tri.size(), 0.01f, idx.data(), &n);
    check(rc == GEOMLIB_MESH_SUCCESS && n == 3, "triangle keeps all corners");

    if (failures != 0) {
        std::printf("%d simplify test(s) failed\n", failures);
        return 1;
    }
    std::printf("All simplify tests passed\n");
    return 0;
}
