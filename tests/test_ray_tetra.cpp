/*
 * Ray-tetrahedron: entry hit inside the solid, exit from an interior ray, and
 * a ray pointing away.
 */
#include <geomlib/geomlib_ray_tetra.h>
#include <geomlib/geomlib_tetramesh.h>

#include <cmath>
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
    geomlib_ray_tetra_hit h;

    geomlib_vec3 o = geomlib_vec3_create(-1, -1, -1);
    geomlib_vec3 d = geomlib_vec3_create(1, 1, 1);
    int hit = geomlib_ray_tetrahedron_intersect(o, d, v0, v1, v2, v3, &h);
    std::printf("ray-tetra hit=%d t=%.4f point=(%.3f,%.3f,%.3f)\n",
                hit, h.t, h.point.x, h.point.y, h.point.z);
    check(hit && h.t > 0.0f &&
          geomlib_point_in_tetrahedron(h.point, v0, v1, v2, v3),
          "ray enters the tetrahedron");

    hit = geomlib_ray_tetrahedron_intersect(
        geomlib_vec3_create(0.25f, 0.25f, 0.25f), d, v0, v1, v2, v3, &h);
    check(hit && h.t > 0.0f, "interior ray exits the tetrahedron");

    hit = geomlib_ray_tetrahedron_intersect(
        geomlib_vec3_create(3, 3, 3), d, v0, v1, v2, v3, &h);
    check(!hit, "ray pointing away misses");
    if (failures != 0) {
        std::printf("%d ray-tetrahedron test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ray-tetrahedron tests passed\n");
    return 0;
}
