/*
 * Ray-cylinder: radial hit, interior origin, miss, and parallel ray.
 */
#include <geomlib/geomlib_ray_cylinder.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 o = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 axis_dir = geomlib_vec3_create(0, 0, 1);
    geomlib_ray_cylinder_hit h;

    /* Unit cylinder on the z axis, ray from (-2,0,0) along +x. */
    int hit = geomlib_ray_cylinder_intersect(
        geomlib_vec3_create(-2, 0, 0), geomlib_vec3_create(1, 0, 0),
        o, axis_dir, 1.0f, &h);
    std::printf("ray-cylinder hit=%d t=%.4f point=(%.2f,%.2f,%.2f) n=(%.2f,%.2f,%.2f)\n",
                hit, h.t, h.point.x, h.point.y, h.point.z,
                h.normal.x, h.normal.y, h.normal.z);
    check(hit && std::fabs(h.t - 1.0) < 1e-6 &&
          std::fabs(h.point.x + 1.0f) < 1e-5f &&
          std::fabs(h.normal.x + 1.0f) < 1e-5f,
          "entry point and outward normal are exact");

    /* Ray starting inside exits at +x. */
    hit = geomlib_ray_cylinder_intersect(
        o, geomlib_vec3_create(1, 0, 0), o, axis_dir, 1.0f, &h);
    check(hit && std::fabs(h.t - 1.0) < 1e-6 &&
          std::fabs(h.normal.x - 1.0f) < 1e-5f,
          "ray from inside exits at the far wall");

    /* Parallel ray offset outside: no hit. */
    hit = geomlib_ray_cylinder_intersect(
        geomlib_vec3_create(2, 0, -2), geomlib_vec3_create(0, 0, 1),
        o, axis_dir, 1.0f, &h);
    check(!hit, "parallel ray outside misses");

    /* Radially offset ray that misses. */
    hit = geomlib_ray_cylinder_intersect(
        geomlib_vec3_create(-2, 2, 0), geomlib_vec3_create(1, 0, 0),
        o, axis_dir, 1.0f, &h);
    check(!hit, "offset ray misses the cylinder");
    if (failures != 0) {
        std::printf("%d ray-cylinder test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ray-cylinder tests passed\n");
    return 0;
}
