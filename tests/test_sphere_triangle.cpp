/*
 * Sphere-triangle: interior overlap, edge-proximity penetration and miss.
 */
#include <geomlib/geomlib_sphere_triangle.h>

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
    geomlib_vec3 v1 = geomlib_vec3_create(2, 0, 0);
    geomlib_vec3 v2 = geomlib_vec3_create(0, 2, 0);
    geomlib_sphere_triangle_result r;

    /* Just above the triangle interior. */
    int hit = geomlib_sphere_triangle_overlap(
        geomlib_vec3_create(0.5f, 0.5f, 0.1f), 0.2f, v0, v1, v2, &r);
    std::printf("sphere-tri hit=%d closest=(%.2f,%.2f,%.2f) d=%.4f\n",
                hit, r.closest.x, r.closest.y, r.closest.z, r.distance);
    check(hit && std::fabs(r.closest.z) < 1e-4f &&
          std::fabs(r.distance + 0.1f) < 1e-4f,
          "sphere penetrates the triangle interior");

    /* Near the hypotenuse edge: closest point lies on x+y=2. */
    hit = geomlib_sphere_triangle_overlap(
        geomlib_vec3_create(1, 1, -2), 2.5f, v0, v1, v2, &r);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-4f &&
          std::fabs(r.closest.x + r.closest.y - 2.0f) < 1e-3f,
          "closest point on the edge gives penetration 0.5");

    hit = geomlib_sphere_triangle_overlap(
        geomlib_vec3_create(5, 5, 5), 0.5f, v0, v1, v2, &r);
    check(!hit, "distant sphere misses the triangle");
    if (failures != 0) {
        std::printf("%d sphere-triangle test(s) failed\n", failures);
        return 1;
    }
    std::printf("All sphere-triangle tests passed\n");
    return 0;
}
