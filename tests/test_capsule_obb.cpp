/*
 * Capsule-OBB: axis-aligned gap and the same gap under rotation.
 */
#include <geomlib/geomlib_capsule_obb.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_obb box;
    geomlib_capsule_obb_result r;
    box.center = geomlib_vec3_create(0, 0, 0);
    box.axes[0] = geomlib_vec3_create(1, 0, 0);
    box.axes[1] = geomlib_vec3_create(0, 1, 0);
    box.axes[2] = geomlib_vec3_create(0, 0, 1);
    box.half_extents = geomlib_vec3_create(1, 1, 1);

    int hit = geomlib_capsule_obb_overlap(
        geomlib_vec3_create(2, 0.5f, 0.5f), geomlib_vec3_create(3, 0.5f, 0.5f),
        0.5f, &box, &r);
    std::printf("capsule-obb miss d=%.4f\n", r.distance);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-4f,
          "axis one unit from the box misses with radius 0.5");

    hit = geomlib_capsule_obb_overlap(
        geomlib_vec3_create(2, 0.5f, 0.5f), geomlib_vec3_create(3, 0.5f, 0.5f),
        1.5f, &box, &r);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-4f,
          "radius larger than the gap penetrates");

    float s = (float)(std::sqrt(2.0) / 2.0);
    box.axes[0] = geomlib_vec3_create(s, s, 0);
    box.axes[1] = geomlib_vec3_create(-s, s, 0);
    box.half_extents = geomlib_vec3_create(1, 0.5f, 0.5f);
    hit = geomlib_capsule_obb_overlap(
        geomlib_vec3_create(2 * s, -2 * s, 0),
        geomlib_vec3_create(1.5f * s, -1.5f * s, 0), 0.5f, &box, &r);
    std::printf("rotated capsule-obb d=%.4f\n", r.distance);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-3f,
          "rotated box uses local coordinates");
    if (failures != 0) {
        std::printf("%d capsule-OBB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All capsule-OBB tests passed\n");
    return 0;
}
