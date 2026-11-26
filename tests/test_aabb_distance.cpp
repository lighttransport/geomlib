/*
 * AABB distance: axis gap, overlap and diagonal separation with closest points.
 */
#include <geomlib/geomlib_aabb_distance.h>

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
    geomlib_vec3 a1 = geomlib_vec3_create(1, 1, 1);
    geomlib_aabb_distance_result r;

    int ok = geomlib_aabb_distance(
        a0, a1, geomlib_vec3_create(2, 0, 0), geomlib_vec3_create(3, 1, 1), &r);
    std::printf("aabb dist=%.4f closest=(%.2f,%.2f,%.2f)-(%.2f,%.2f,%.2f)\n",
                r.distance, r.point_a.x, r.point_a.y, r.point_a.z,
                r.point_b.x, r.point_b.y, r.point_b.z);
    check(ok && std::fabs(r.distance - 1.0f) < 1e-5f &&
          std::fabs(r.point_a.x - 1.0f) < 1e-5f &&
          std::fabs(r.point_b.x - 2.0f) < 1e-5f,
          "unit gap along x with correct closest points");

    ok = geomlib_aabb_distance(
        a0, a1, geomlib_vec3_create(0.5f, 0.5f, 0.5f),
        geomlib_vec3_create(2, 2, 2), &r);
    check(ok && r.distance == 0.0f, "overlapping boxes have distance zero");

    ok = geomlib_aabb_distance(
        a0, a1, geomlib_vec3_create(2, 2, 0), geomlib_vec3_create(3, 3, 1), &r);
    std::printf("aabb diagonal dist=%.6f expect=%.6f\n", r.distance, std::sqrt(2.0));
    check(ok && std::fabs(r.distance - std::sqrt(2.0)) < 1e-4,
          "diagonal gap combines axes");
    if (failures != 0) {
        std::printf("%d AABB distance test(s) failed\n", failures);
        return 1;
    }
    std::printf("All AABB distance tests passed\n");
    return 0;
}
