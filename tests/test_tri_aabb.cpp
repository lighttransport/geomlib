/*
 * Triangle-AABB: crossing, contained and separated cases.
 */
#include <geomlib/geomlib_tri_aabb.h>

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

    /* Triangle crossing the unit box. */
    int hit = geomlib_triangle_aabb_overlap(
        geomlib_vec3_create(0.5f, 0.5f, -1),
        geomlib_vec3_create(0.5f, 1.5f, 0.5f),
        geomlib_vec3_create(1.5f, 0.5f, 0.5f), lo, hi);
    std::printf("tri-aabb crossing hit=%d\n", hit);
    check(hit == 1, "crossing triangle overlaps the box");

    /* Triangle entirely inside. */
    hit = geomlib_triangle_aabb_overlap(
        geomlib_vec3_create(0.2f, 0.2f, 0.2f),
        geomlib_vec3_create(0.8f, 0.2f, 0.2f),
        geomlib_vec3_create(0.2f, 0.8f, 0.2f), lo, hi);
    check(hit == 1, "interior triangle overlaps the box");

    /* Triangle above the box. */
    hit = geomlib_triangle_aabb_overlap(
        geomlib_vec3_create(0.2f, 0.2f, 3),
        geomlib_vec3_create(0.8f, 0.2f, 3),
        geomlib_vec3_create(0.2f, 0.8f, 3), lo, hi);
    check(hit == 0, "triangle above the box is separated");

    /* Triangle beside the box but its projected area overlaps: still outside. */
    hit = geomlib_triangle_aabb_overlap(
        geomlib_vec3_create(2, 0.2f, 0.2f),
        geomlib_vec3_create(2, 0.8f, 0.2f),
        geomlib_vec3_create(2, 0.2f, 0.8f), lo, hi);
    check(hit == 0, "triangle beside the box is separated");
    if (failures != 0) {
        std::printf("%d tri-AABB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All tri-AABB tests passed\n");
    return 0;
}
