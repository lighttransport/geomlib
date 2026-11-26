/*
 * Triangle-triangle intersection: coplanar crossing segment, disjoint case
 * and perpendicular case.
 */
#include <geomlib/geomlib_tri_tri.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 a0 = geomlib_vec3_create(-1, -1, 0);
    geomlib_vec3 a1 = geomlib_vec3_create(1, -1, 0);
    geomlib_vec3 a2 = geomlib_vec3_create(0, 1, 0);
    geomlib_vec3 b0 = geomlib_vec3_create(-1, 1, 0);
    geomlib_vec3 b1 = geomlib_vec3_create(1, 1, 0);
    geomlib_vec3 b2 = geomlib_vec3_create(0, -1, 0);

    geomlib_tri_tri_result r = geomlib_triangle_intersection(a0, a1, a2, b0, b1, b2);
    std::printf("tri tri hit=%d seg=(%.3f,%.3f)-(%.3f,%.3f)\n",
                r.intersects, r.p0.x, r.p0.y, r.p1.x, r.p1.y);
    float ylo = r.p0.y < r.p1.y ? r.p0.y : r.p1.y;
    float yhi = r.p0.y < r.p1.y ? r.p1.y : r.p0.y;
    check(r.intersects == 1 &&
          std::fabs(r.p0.x) < 1e-3f && std::fabs(r.p1.x) < 1e-3f &&
          std::fabs(ylo + 1.0f) < 1e-3f && std::fabs(yhi - 1.0f) < 1e-3f,
          "coplanar crossing triangles intersect");
    check(std::fabs(r.p0.z) < 1e-4f && std::fabs(r.p1.z) < 1e-4f,
          "segment lies in the plane");

    geomlib_vec3 up = geomlib_vec3_create(0, 0, 5);
    r = geomlib_triangle_intersection(
        a0, a1, a2,
        geomlib_vec3_add(b0, up), geomlib_vec3_add(b1, up),
        geomlib_vec3_add(b2, up));
    check(r.intersects == 0, "disjoint triangles do not intersect");

    r = geomlib_triangle_intersection(
        a0, a1, a2,
        geomlib_vec3_create(0, -2, -2), geomlib_vec3_create(0, -2, 2),
        geomlib_vec3_create(0, 2, 0));
    std::printf("tri tri perp hit=%d\n", r.intersects);
    check(r.intersects == 1, "perpendicular triangles intersect");

    if (failures != 0) {
        std::printf("%d tri-tri test(s) failed\n", failures);
        return 1;
    }
    std::printf("All tri-tri tests passed\n");
    return 0;
}
