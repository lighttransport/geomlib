/*
 * Capsule: tangential, penetrating, separated and crossing capsules; signed
 * point-to-capsule distance.
 */
#include <geomlib/geomlib_capsule.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_capsule_result r;

    /* Collinear capsules 1 apart with radius 0.5 each: tangent. */
    int hit = geomlib_capsule_capsule_overlap(
        geomlib_vec3_create(0, 0, 0), geomlib_vec3_create(0, 0, 2), 0.5f,
        geomlib_vec3_create(0, 0, 3), geomlib_vec3_create(0, 0, 5), 0.5f, &r);
    std::printf("capsule tangent hit=%d d=%.4f\n", hit, r.distance);
    check(hit && std::fabs(r.distance) < 1e-4f,
          "collinear capsules touching are reported");

    /* Move the second capsule away: gap 2. */
    hit = geomlib_capsule_capsule_overlap(
        geomlib_vec3_create(0, 0, 0), geomlib_vec3_create(0, 0, 2), 0.5f,
        geomlib_vec3_create(0, 0, 5), geomlib_vec3_create(0, 0, 7), 0.5f, &r);
    check(!hit && std::fabs(r.distance - 2.0f) < 1e-4f,
          "separated capsules report the gap");

    /* Perpendicular axes crossing at (1,0,0): penetration depth 1. */
    hit = geomlib_capsule_capsule_overlap(
        geomlib_vec3_create(0, 0, 0), geomlib_vec3_create(2, 0, 0), 0.5f,
        geomlib_vec3_create(1, -1, 0), geomlib_vec3_create(1, 1, 0), 0.5f, &r);
    check(hit && std::fabs(r.distance + 1.0f) < 1e-4f,
          "crossing capsules penetrate by the sum of radii");

    /* Point distance to a z-axis capsule of radius 0.5. */
    float d = geomlib_capsule_point_distance(
        geomlib_vec3_create(2, 0, 0), geomlib_vec3_create(0, 0, 0),
        geomlib_vec3_create(0, 0, 2), 0.5f, nullptr);
    std::printf("capsule point distance=%.4f expect=1.5\n", d);
    check(std::fabs(d - 1.5f) < 1e-4f, "point distance is exact");
    if (failures != 0) {
        std::printf("%d capsule test(s) failed\n", failures);
        return 1;
    }
    std::printf("All capsule tests passed\n");
    return 0;
}
