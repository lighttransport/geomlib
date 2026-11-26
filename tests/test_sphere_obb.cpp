/*
 * Sphere-OBB: axis-aligned miss/hit and a rotated box case.
 */
#include <geomlib/geomlib_sphere_obb.h>

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
    geomlib_sphere_obb_result r;
    box.center = geomlib_vec3_create(0, 0, 0);
    box.axes[0] = geomlib_vec3_create(1, 0, 0);
    box.axes[1] = geomlib_vec3_create(0, 1, 0);
    box.axes[2] = geomlib_vec3_create(0, 0, 1);
    box.half_extents = geomlib_vec3_create(1, 1, 1);

    int hit = geomlib_sphere_obb_overlap(
        geomlib_vec3_create(2, 0, 0), 0.5f, &box, &r);
    std::printf("sphere-obb miss d=%.4f closest=(%.2f,%.2f,%.2f)\n",
                r.distance, r.closest.x, r.closest.y, r.closest.z);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-5f &&
          std::fabs(r.closest.x - 1.0f) < 1e-5f,
          "sphere outside the box misses with face projection");

    hit = geomlib_sphere_obb_overlap(
        geomlib_vec3_create(2, 0, 0), 1.5f, &box, &r);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-5f,
          "sphere reaching the box penetrates");

    float s = (float)(std::sqrt(2.0) / 2.0);
    box.axes[0] = geomlib_vec3_create(s, s, 0);
    box.axes[1] = geomlib_vec3_create(-s, s, 0);
    box.half_extents = geomlib_vec3_create(1, 0.5f, 0.5f);
    hit = geomlib_sphere_obb_overlap(
        geomlib_vec3_create(2 * s, 2 * s, 0), 0.5f, &box, &r);
    std::printf("rotated sphere-obb d=%.4f\n", r.distance);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-4f,
          "rotated box uses local coordinates");
    if (failures != 0) {
        std::printf("%d sphere-OBB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All sphere-OBB tests passed\n");
    return 0;
}
