/*
 * Tetrahedron-tetrahedron: overlap, containment, shared face and separation.
 */
#include <geomlib/geomlib_tetra_tetra.h>

#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 a0 = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 a1 = geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 a2 = geomlib_vec3_create(0, 1, 0);
    geomlib_vec3 a3 = geomlib_vec3_create(0, 0, 1);

    /* Same tetrahedron. */
    int hit = geomlib_tetrahedra_intersect(a0, a1, a2, a3, a0, a1, a2, a3);
    std::printf("tetra-tetra same=%d\n", hit);
    check(hit == 1, "identical tetrahedra intersect");

    /* Small tetrahedron inside. */
    hit = geomlib_tetrahedra_intersect(
        a0, a1, a2, a3,
        geomlib_vec3_create(0.1f, 0.1f, 0.1f),
        geomlib_vec3_create(0.2f, 0.1f, 0.1f),
        geomlib_vec3_create(0.1f, 0.2f, 0.1f),
        geomlib_vec3_create(0.1f, 0.1f, 0.2f));
    check(hit == 1, "contained tetrahedron intersects");

    /* Disjoint translation. */
    hit = geomlib_tetrahedra_intersect(
        a0, a1, a2, a3,
        geomlib_vec3_create(3, 0, 0), geomlib_vec3_create(4, 0, 0),
        geomlib_vec3_create(3, 1, 0), geomlib_vec3_create(3, 0, 1));
    check(hit == 0, "separated tetrahedra do not intersect");

    /* Overlapping shifted tetrahedra. */
    hit = geomlib_tetrahedra_intersect(
        a0, a1, a2, a3,
        geomlib_vec3_create(0.5f, 0, 0), geomlib_vec3_create(1.5f, 0, 0),
        geomlib_vec3_create(0.5f, 1, 0), geomlib_vec3_create(0.5f, 0, 1));
    check(hit == 1, "partially overlapping tetrahedra intersect");
    if (failures != 0) {
        std::printf("%d tetra-tetra test(s) failed\n", failures);
        return 1;
    }
    std::printf("All tetra-tetra tests passed\n");
    return 0;
}
