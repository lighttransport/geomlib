/*
 * Plane-sphere: great circle, offset circle, tangency and disjoint case.
 */
#include <geomlib/geomlib_plane_sphere.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 p = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 n = geomlib_vec3_create(0, 0, 1);

    geomlib_plane_sphere_intersection r = geomlib_plane_sphere_intersect(
        p, n, geomlib_vec3_create(0, 0, 0), 2.0f);
    std::printf("plane-sphere great circle kind=%d radius=%.4f\n",
                (int)r.kind, r.radius);
    check(r.kind == GEOMLIB_PLANE_SPHERE_CIRCLE &&
          std::fabs(r.radius - 2.0f) < 1e-5f &&
          std::fabs(r.center.z) < 1e-5f, "great circle has radius 2");

    r = geomlib_plane_sphere_intersect(
        p, n, geomlib_vec3_create(0, 0, 1), 2.0f);
    std::printf("offset circle radius=%.6f expect=%.6f\n",
                r.radius, std::sqrt(3.0));
    check(r.kind == GEOMLIB_PLANE_SPHERE_CIRCLE &&
          std::fabs(r.radius - std::sqrt(3.0)) < 1e-5f &&
          std::fabs(r.center.z) < 1e-5f, "offset circle is at the origin");

    r = geomlib_plane_sphere_intersect(
        p, n, geomlib_vec3_create(0, 0, 2), 2.0f);
    check(r.kind == GEOMLIB_PLANE_SPHERE_TANGENT, "tangent sphere");

    r = geomlib_plane_sphere_intersect(
        p, n, geomlib_vec3_create(0, 0, 3), 2.0f);
    check(r.kind == GEOMLIB_PLANE_SPHERE_NONE, "disjoint sphere");
    if (failures != 0) {
        std::printf("%d plane-sphere test(s) failed\n", failures);
        return 1;
    }
    std::printf("All plane-sphere tests passed\n");
    return 0;
}
