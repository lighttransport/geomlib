/*
 * Triangle-triangle distance: separated, intersecting and parallel-plane cases.
 */
#include <geomlib/geomlib_tri_tri_distance.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 a0 = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 a1 = geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 a2 = geomlib_vec3_create(0, 1, 0);
    geomlib_tri_tri_distance_result r;

    /* Shifted 3 along x: closest vertices (1,0,0) and (3,0,0). */
    int ok = geomlib_triangle_triangle_distance(
        a0, a1, a2,
        geomlib_vec3_create(3, 0, 0), geomlib_vec3_create(4, 0, 0),
        geomlib_vec3_create(3, 1, 0), &r);
    std::printf("tri-tri dist=%.6f\n", r.distance);
    check(ok && std::fabs(r.distance - 2.0f) < 1e-4f &&
          std::fabs(r.point_a.x - 1.0f) < 1e-4f &&
          std::fabs(r.point_b.x - 3.0f) < 1e-4f,
          "separated triangles have vertex-to-vertex distance 2");

    /* Overlapping coplanar triangles. */
    ok = geomlib_triangle_triangle_distance(
        a0, a1, a2,
        geomlib_vec3_create(0.5f, 0.5f, 0), geomlib_vec3_create(0.8f, 0.5f, 0),
        geomlib_vec3_create(0.5f, 0.8f, 0), &r);
    check(ok && r.distance == 0.0f, "intersecting triangles have distance zero");

    /* Parallel planes two apart with overlapping projection. */
    ok = geomlib_triangle_triangle_distance(
        a0, a1, a2,
        geomlib_vec3_create(0, 0, 2), geomlib_vec3_create(1, 0, 2),
        geomlib_vec3_create(0, 1, 2), &r);
    std::printf("parallel-plane dist=%.6f\n", r.distance);
    check(ok && std::fabs(r.distance - 2.0f) < 1e-4f,
          "parallel triangles are two apart");
    if (failures != 0) {
        std::printf("%d tri-tri distance test(s) failed\n", failures);
        return 1;
    }
    std::printf("All tri-tri distance tests passed\n");
    return 0;
}
