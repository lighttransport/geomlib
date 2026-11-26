/*
 * Sphere-capsule: tangency, miss, and penetration at the axis midpoint.
 */
#include <geomlib/geomlib_sphere_capsule.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 a = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 b = geomlib_vec3_create(2, 0, 0);
    geomlib_sphere_capsule_result r;

    /* Axis one unit from the sphere centre, radii summing to one: tangent. */
    int hit = geomlib_sphere_capsule_overlap(
        geomlib_vec3_create(1, 0, 1), 0.5f, a, b, 0.5f, &r);
    std::printf("sphere-capsule tangent d=%.4f point=(%.2f,%.2f,%.2f)\n",
                r.distance, r.point.x, r.point.y, r.point.z);
    check(hit && std::fabs(r.distance) < 1e-5f &&
          std::fabs(r.point.x - 1.0f) < 1e-5f,
          "radii summing to the gap are tangent at the closest axis point");

    hit = geomlib_sphere_capsule_overlap(
        geomlib_vec3_create(1, 0, 1), 0.2f, a, b, 0.5f, &r);
    check(!hit && std::fabs(r.distance - 0.3f) < 1e-5f,
          "smaller sphere misses with gap 0.3");

    hit = geomlib_sphere_capsule_overlap(
        geomlib_vec3_create(1, 0, 0), 0.5f, a, b, 0.5f, &r);
    check(hit && std::fabs(r.distance + 1.0f) < 1e-5f,
          "centre on the axis penetrates by the sum of radii");
    if (failures != 0) {
        std::printf("%d sphere-capsule test(s) failed\n", failures);
        return 1;
    }
    std::printf("All sphere-capsule tests passed\n");
    return 0;
}
