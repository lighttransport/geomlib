/*
 * Triangle-tetrahedron: crossing, inside, vertex-in-triangle and separated.
 */
#include <geomlib/geomlib_tri_tetra.h>

#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 v0 = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 v1 = geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 v2 = geomlib_vec3_create(0, 1, 0);
    geomlib_vec3 v3 = geomlib_vec3_create(0, 0, 1);

    /* Triangle cutting through the tetrahedron. */
    int hit = geomlib_triangle_tetrahedron_intersect(
        geomlib_vec3_create(-1, 0.3f, 0.3f), geomlib_vec3_create(2, 0.3f, 0.3f),
        geomlib_vec3_create(-1, 0.8f, 0.3f), v0, v1, v2, v3);
    std::printf("tri-tetra crossing=%d\n", hit);
    check(hit == 1, "cutting triangle intersects");

    /* Triangle entirely inside. */
    hit = geomlib_triangle_tetrahedron_intersect(
        geomlib_vec3_create(0.1f, 0.1f, 0.1f),
        geomlib_vec3_create(0.3f, 0.1f, 0.1f),
        geomlib_vec3_create(0.1f, 0.3f, 0.1f), v0, v1, v2, v3);
    check(hit == 1, "interior triangle intersects");

    /* Tetrahedron vertex inside an external triangle: use the face x+y+z=1
     * triangle itself, which shares vertices with the solid. */
    hit = geomlib_triangle_tetrahedron_intersect(v1, v2, v3, v0, v1, v2, v3);
    check(hit == 1, "shared face intersects");

    /* Far away. */
    hit = geomlib_triangle_tetrahedron_intersect(
        geomlib_vec3_create(5, 5, 5), geomlib_vec3_create(6, 5, 5),
        geomlib_vec3_create(5, 6, 5), v0, v1, v2, v3);
    check(hit == 0, "distant triangle misses");
    if (failures != 0) {
        std::printf("%d tri-tetra test(s) failed\n", failures);
        return 1;
    }
    std::printf("All tri-tetra tests passed\n");
    return 0;
}
