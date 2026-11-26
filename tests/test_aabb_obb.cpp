/*
 * AABB-OBB: identical boxes, separation and a rotated overlap.
 */
#include <geomlib/geomlib_aabb_obb.h>

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
    geomlib_vec3 hi = geomlib_vec3_create(2, 2, 2);
    geomlib_obb box;
    box.center = geomlib_vec3_create(1, 1, 1);
    box.axes[0] = geomlib_vec3_create(1, 0, 0);
    box.axes[1] = geomlib_vec3_create(0, 1, 0);
    box.axes[2] = geomlib_vec3_create(0, 0, 1);
    box.half_extents = geomlib_vec3_create(1, 1, 1);

    int hit = geomlib_aabb_obb_overlap(lo, hi, &box);
    std::printf("aabb-obb identical=%d\n", hit);
    check(hit == 1, "coincident AABB and OBB overlap");

    box.center = geomlib_vec3_create(6, 1, 1);
    hit = geomlib_aabb_obb_overlap(lo, hi, &box);
    check(hit == 0, "separated boxes do not overlap");

    float s = (float)(std::sqrt(2.0) / 2.0);
    box.center = geomlib_vec3_create(1, 1, 1);
    box.axes[0] = geomlib_vec3_create(s, s, 0);
    box.axes[1] = geomlib_vec3_create(-s, s, 0);
    box.half_extents = geomlib_vec3_create(1, 0.5f, 1);
    hit = geomlib_aabb_obb_overlap(lo, hi, &box);
    check(hit == 1, "rotated overlapping box is detected");
    if (failures != 0) {
        std::printf("%d AABB-OBB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All AABB-OBB tests passed\n");
    return 0;
}
