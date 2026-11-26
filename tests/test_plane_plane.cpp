/*
 * Plane-plane: perpendicular planes meet in a line through the right point,
 * points on the line satisfy both plane equations, parallel planes report no
 * intersection.
 */
#include <geomlib/geomlib_plane_plane.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 p0 = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 n0 = geomlib_vec3_create(0, 0, 1);
    geomlib_vec3 p1 = geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 n1 = geomlib_vec3_create(1, 0, 0);
    geomlib_plane_line l;

    int hit = geomlib_plane_plane_intersection(p0, n0, p1, n1, &l);
    std::printf("plane-plane hit=%d point=(%.2f,%.2f,%.2f) dir=(%.2f,%.2f,%.2f)\n",
                hit, l.point.x, l.point.y, l.point.z,
                l.direction.x, l.direction.y, l.direction.z);
    check(hit && std::fabs(l.point.x - 1.0f) < 1e-4f &&
          std::fabs(l.point.z) < 1e-4f, "line passes through x=1, z=0");
    check(std::fabs(l.direction.x) < 1e-4f &&
          std::fabs(l.direction.y - 1.0f) < 1e-4f &&
          std::fabs(l.direction.z) < 1e-4f, "direction is the y axis");

    /* Points on the line satisfy both plane equations. */
    geomlib_vec3 q = geomlib_vec3_add(l.point, geomlib_vec3_mul(l.direction, 3.0f));
    check(std::fabs(q.x - 1.0f) < 1e-3f && std::fabs(q.z) < 1e-3f,
          "points on the line satisfy both planes");

    /* Parallel planes. */
    hit = geomlib_plane_plane_intersection(
        p0, n0, geomlib_vec3_create(0, 0, 5), n0, &l);
    check(!hit && l.parallel, "parallel planes report no line");
    if (failures != 0) {
        std::printf("%d plane-plane test(s) failed\n", failures);
        return 1;
    }
    std::printf("All plane-plane tests passed\n");
    return 0;
}
