/*
 * OBB-OBB: identical, separated and rotated overlap.
 */
#include <geomlib/geomlib_obb_obb.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static void axis_aligned(geomlib_obb *b, geomlib_vec3 c)
{
    b->center = c;
    b->axes[0] = geomlib_vec3_create(1, 0, 0);
    b->axes[1] = geomlib_vec3_create(0, 1, 0);
    b->axes[2] = geomlib_vec3_create(0, 0, 1);
    b->half_extents = geomlib_vec3_create(1, 1, 1);
}

int main(void)
{
    geomlib_obb a, b;
    axis_aligned(&a, geomlib_vec3_create(0, 0, 0));
    axis_aligned(&b, geomlib_vec3_create(0, 0, 0));

    int hit = geomlib_obb_obb_overlap(&a, &b);
    std::printf("obb-obb identical=%d\n", hit);
    check(hit == 1, "coincident OBBs overlap");

    axis_aligned(&b, geomlib_vec3_create(5, 0, 0));
    check(geomlib_obb_obb_overlap(&a, &b) == 0, "separated OBBs do not overlap");

    float s = (float)(std::sqrt(2.0) / 2.0);
    axis_aligned(&b, geomlib_vec3_create(0, 0, 0));
    b.axes[0] = geomlib_vec3_create(s, s, 0);
    b.axes[1] = geomlib_vec3_create(-s, s, 0);
    hit = geomlib_obb_obb_overlap(&a, &b);
    check(hit == 1, "rotated overlapping OBB is detected");
    if (failures != 0) {
        std::printf("%d OBB-OBB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All OBB-OBB tests passed\n");
    return 0;
}
