/*
 * Arc length: cumulative table, total length and interpolation.
 */
#include <geomlib/geomlib_arc_length.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 p[3] = {
        geomlib_vec3_create(0, 0, 0),
        geomlib_vec3_create(3, 0, 0),
        geomlib_vec3_create(3, 4, 0)
    };
    double table[3], total;
    int rc = geomlib_polyline_arc_table(p, 3, table, &total);
    std::printf("arc table rc=%d cum=(%.2f,%.2f,%.2f) total=%.2f\n",
                rc, table[0], table[1], table[2], total);
    check(rc == GEOMLIB_MESH_SUCCESS && std::fabs(total - 7.0) < 1e-9 &&
          std::fabs(table[1] - 3.0) < 1e-9,
          "cumulative lengths and total are exact");

    geomlib_vec3 q = geomlib_polyline_at_length(p, 3, table, 1.5);
    check(std::fabs(q.x - 1.5f) < 1e-5f && std::fabs(q.y) < 1e-5f,
          "interpolation inside the first segment");
    q = geomlib_polyline_at_length(p, 3, table, 5.0);
    check(std::fabs(q.x - 3.0f) < 1e-5f && std::fabs(q.y - 2.0f) < 1e-5f,
          "interpolation inside the second segment");
    q = geomlib_polyline_at_length(p, 3, table, 100.0);
    check(std::fabs(q.y - 4.0f) < 1e-5f, "over-range clamps to the end");
    check(std::fabs(geomlib_polyline_length(p, 3) - 7.0) < 1e-9,
          "total length helper");
    if (failures != 0) {
        std::printf("%d arc length test(s) failed\n", failures);
        return 1;
    }
    std::printf("All arc length tests passed\n");
    return 0;
}
