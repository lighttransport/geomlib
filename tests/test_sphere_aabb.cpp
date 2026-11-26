/*
 * Sphere/AABB: tangency, penetration, interior centre and a clean miss.
 */
#include <geomlib/geomlib_sphere_aabb.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 lo = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 hi = geomlib_vec3_create(1, 1, 1);
    geomlib_sphere_aabb_result r;

    /* Unit distance to the box face: radius 0.5 misses, 1.5 penetrates. */
    geomlib_vec3 c = geomlib_vec3_create(2, 0.5f, 0.5f);
    int hit = geomlib_sphere_aabb_overlap(c, 0.5f, lo, hi, &r);
    std::printf("sphere-aabb miss d=%.4f\n", r.distance);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-5f,
          "sphere outside does not overlap");
    check(std::fabs(r.closest.x - 1.0f) < 1e-5f &&
          std::fabs(r.closest.y - 0.5f) < 1e-5f,
          "closest box point is the face projection");

    hit = geomlib_sphere_aabb_overlap(c, 1.5f, lo, hi, &r);
    std::printf("sphere-aabb hit d=%.4f\n", r.distance);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-5f,
          "sphere that reaches the box overlaps");

    /* Interior centre. */
    hit = geomlib_sphere_aabb_overlap(
        geomlib_vec3_create(0.5f, 0.5f, 0.5f), 0.2f, lo, hi, &r);
    check(hit && r.distance < 0.0f, "interior centre overlaps");
    if (failures != 0) {
        std::printf("%d sphere-AABB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All sphere-AABB tests passed\n");
    return 0;
}
