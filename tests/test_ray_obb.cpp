/*
 * Ray-OBB: axis-aligned hit, rotated hit, interior origin and a miss.
 */
#include <geomlib/geomlib_ray_obb.h>

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
    geomlib_ray_obb_hit h;
    box.center = geomlib_vec3_create(0, 0, 0);
    box.axes[0] = geomlib_vec3_create(1, 0, 0);
    box.axes[1] = geomlib_vec3_create(0, 1, 0);
    box.axes[2] = geomlib_vec3_create(0, 0, 1);
    box.half_extents = geomlib_vec3_create(1, 1, 1);

    int hit = geomlib_ray_obb_intersect(
        geomlib_vec3_create(-2, 0, 0), geomlib_vec3_create(1, 0, 0), &box, &h);
    std::printf("ray-obb hit=%d t=%.4f point=(%.2f,%.2f,%.2f) n=(%.2f,%.2f,%.2f)\n",
                hit, h.t, h.point.x, h.point.y, h.point.z,
                h.normal.x, h.normal.y, h.normal.z);
    check(hit && std::fabs(h.t - 1.0) < 1e-6 &&
          std::fabs(h.point.x + 1.0f) < 1e-5f &&
          std::fabs(h.normal.x + 1.0f) < 1e-5f,
          "axis-aligned entry and outward normal");

    /* Rotated box: local y face at -0.5 along (-0.707,0.707,0). */
    float s = (float)(std::sqrt(2.0) / 2.0);
    box.axes[0] = geomlib_vec3_create(s, s, 0);
    box.axes[1] = geomlib_vec3_create(-s, s, 0);
    box.half_extents = geomlib_vec3_create(1, 0.5f, 0.5f);
    hit = geomlib_ray_obb_intersect(
        geomlib_vec3_create(0, -2, 0), geomlib_vec3_create(0, 1, 0), &box, &h);
    check(hit && std::fabs((h.point.x - box.center.x) * box.axes[1].x +
                           (h.point.y - box.center.y) * box.axes[1].y + 0.5f) < 1e-4f,
          "rotated box hit lies on the correct local face");

    /* Ray from inside exits through +y local face. */
    hit = geomlib_ray_obb_intersect(
        box.center, geomlib_vec3_create(0, 1, 0), &box, &h);
    check(hit && h.t > 0.0 && h.normal.y > 0.0f,
          "ray from inside exits with outward normal");

    hit = geomlib_ray_obb_intersect(
        geomlib_vec3_create(5, 0, 0), geomlib_vec3_create(1, 0, 0), &box, &h);
    check(!hit, "ray pointing away misses");
    if (failures != 0) {
        std::printf("%d ray-OBB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ray-OBB tests passed\n");
    return 0;
}
