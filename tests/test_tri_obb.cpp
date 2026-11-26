/*
 * Triangle-OBB: crossing, separated and rotated-box cases.
 */
#include <geomlib/geomlib_tri_obb.h>

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
    box.center = geomlib_vec3_create(0, 0, 0);
    box.axes[0] = geomlib_vec3_create(1, 0, 0);
    box.axes[1] = geomlib_vec3_create(0, 1, 0);
    box.axes[2] = geomlib_vec3_create(0, 0, 1);
    box.half_extents = geomlib_vec3_create(1, 1, 1);

    int hit = geomlib_triangle_obb_overlap(
        geomlib_vec3_create(0.5f, 0.5f, -2),
        geomlib_vec3_create(0.5f, 1.5f, 0.5f),
        geomlib_vec3_create(1.5f, 0.5f, 0.5f), &box);
    std::printf("tri-obb crossing=%d\n", hit);
    check(hit == 1, "triangle crossing the box overlaps");

    hit = geomlib_triangle_obb_overlap(
        geomlib_vec3_create(3, 0.2f, 0.2f),
        geomlib_vec3_create(3, 0.8f, 0.2f),
        geomlib_vec3_create(3, 0.2f, 0.8f), &box);
    check(hit == 0, "triangle beside the box is separated");

    float s = (float)(std::sqrt(2.0) / 2.0);
    box.axes[0] = geomlib_vec3_create(s, s, 0);
    box.axes[1] = geomlib_vec3_create(-s, s, 0);
    box.half_extents = geomlib_vec3_create(1, 0.5f, 0.5f);
    hit = geomlib_triangle_obb_overlap(
        geomlib_vec3_create(0, -1.5f, 0),
        geomlib_vec3_create(0.6f, -0.6f, 0),
        geomlib_vec3_create(-0.6f, -0.6f, 0), &box);
    check(hit == 1, "rotated box is handled in local coordinates");
    if (failures != 0) {
        std::printf("%d tri-OBB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All tri-OBB tests passed\n");
    return 0;
}
