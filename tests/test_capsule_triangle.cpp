/*
 * Capsule-triangle: crossing axis, parallel gap, and an edge-aligned axis.
 */
#include <geomlib/geomlib_capsule_triangle.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 v0 = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 v1 = geomlib_vec3_create(2, 0, 0);
    geomlib_vec3 v2 = geomlib_vec3_create(0, 2, 0);
    geomlib_capsule_triangle_result r;

    int hit = geomlib_capsule_triangle_overlap(
        geomlib_vec3_create(0.5f, 0.5f, -2), geomlib_vec3_create(0.5f, 0.5f, 2),
        0.2f, v0, v1, v2, &r);
    std::printf("capsule-tri crossing d=%.4f\n", r.distance);
    check(hit && std::fabs(r.distance + 0.2f) < 1e-4f,
          "crossing axis penetrates by the radius");

    hit = geomlib_capsule_triangle_overlap(
        geomlib_vec3_create(0.5f, 0.5f, 2), geomlib_vec3_create(0.5f, 0.5f, 4),
        0.5f, v0, v1, v2, &r);
    std::printf("capsule-tri parallel d=%.4f\n", r.distance);
    check(!hit && std::fabs(r.distance - 1.5f) < 1e-4f,
          "parallel axis above the triangle misses with gap 1.5");

    hit = geomlib_capsule_triangle_overlap(
        geomlib_vec3_create(0, 0, -1), geomlib_vec3_create(2, 0, 1),
        0.1f, v0, v1, v2, &r);
    check(hit && std::fabs(r.distance + 0.1f) < 1e-4f,
          "axis crossing the triangle edge overlaps");
    if (failures != 0) {
        std::printf("%d capsule-triangle test(s) failed\n", failures);
        return 1;
    }
    std::printf("All capsule-triangle tests passed\n");
    return 0;
}
