/*
 * Ray-cone: side hit at the expected radius, axial ray misses, and a ray
 * parallel to the axis at a large offset hits further up the cone.
 */
#include <geomlib/geomlib_ray_cone.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 apex = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 axis = geomlib_vec3_create(0, 0, 1);
    float a = (float)(std::acos(std::sqrt(2.0) / 2.0)); /* 45 degrees */
    geomlib_ray_cone_hit h;

    /* At z=2 the cone radius is 2; ray from x=-3 hits at x=-2 (t=1). */
    int hit = geomlib_ray_cone_intersect(
        geomlib_vec3_create(-3, 0, 2), geomlib_vec3_create(1, 0, 0),
        apex, axis, a, &h);
    std::printf("ray-cone hit=%d t=%.4f point=(%.3f,%.3f,%.3f)\n",
                hit, h.t, h.point.x, h.point.y, h.point.z);
    check(hit && std::fabs(h.t - 1.0) < 1e-5 &&
          std::fabs(h.point.x + 2.0f) < 1e-4f,
          "45-degree cone hit at x=-2, z=2");

    /* Ray along the axis away from the apex never reaches the surface. */
    hit = geomlib_ray_cone_intersect(
        geomlib_vec3_create(0, 0, 1), geomlib_vec3_create(0, 0, 1),
        apex, axis, a, &h);
    check(!hit, "axial ray misses the cone surface");

    /* Parallel ray at x=3 hits where the cone radius reaches 3, i.e. z=3. */
    hit = geomlib_ray_cone_intersect(
        geomlib_vec3_create(3, 0, 0), geomlib_vec3_create(0, 0, 1),
        apex, axis, a, &h);
    std::printf("parallel ray-cone t=%.4f\n", h.t);
    check(hit && std::fabs(h.t - 3.0) < 1e-4,
          "parallel ray hits at z=3");
    if (failures != 0) {
        std::printf("%d ray-cone test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ray-cone tests passed\n");
    return 0;
}
