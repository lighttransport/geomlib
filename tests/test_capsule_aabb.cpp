/*
 * Capsule-AABB: crossing axis, separated axis with exact gap, and a miss.
 */
#include <geomlib/geomlib_capsule_aabb.h>

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
    geomlib_capsule_aabb_result r;

    /* Axis passes through the box. */
    int hit = geomlib_capsule_aabb_overlap(
        geomlib_vec3_create(0.5f, 0.5f, -2), geomlib_vec3_create(0.5f, 0.5f, 2),
        0.2f, lo, hi, &r);
    std::printf("capsule-aabb crossing d=%.4f\n", r.distance);
    check(hit && std::fabs(r.distance + 0.2f) < 1e-4f,
          "crossing axis overlaps with penetration 0.2");

    /* Axis parallel to x at distance 1 from the box face. */
    float d = geomlib_segment_aabb_distance(
        geomlib_vec3_create(2, 0.5f, 0.5f), geomlib_vec3_create(3, 0.5f, 0.5f),
        lo, hi);
    std::printf("segment-aabb gap=%.4f\n", d);
    check(std::fabs(d - 1.0f) < 1e-4f, "segment-box distance is 1");

    hit = geomlib_capsule_aabb_overlap(
        geomlib_vec3_create(2, 0.5f, 0.5f), geomlib_vec3_create(3, 0.5f, 0.5f),
        0.5f, lo, hi, &r);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-4f,
          "radius smaller than the gap misses");

    hit = geomlib_capsule_aabb_overlap(
        geomlib_vec3_create(2, 0.5f, 0.5f), geomlib_vec3_create(3, 0.5f, 0.5f),
        1.5f, lo, hi, &r);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-4f,
          "radius larger than the gap penetrates");
    if (failures != 0) {
        std::printf("%d capsule-AABB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All capsule-AABB tests passed\n");
    return 0;
}
