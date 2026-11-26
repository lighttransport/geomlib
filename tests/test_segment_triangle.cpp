/*
 * Segment-triangle: a crossing segment hits at t=0.5 with the right barycentric
 * coordinates; parallel, short and outside triangles miss.
 */
#include <geomlib/geomlib_segment_triangle.h>

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
    geomlib_segment_triangle_hit h;

    int hit = geomlib_segment_triangle_intersect(
        geomlib_vec3_create(0.2f, 0.2f, -1), geomlib_vec3_create(0.2f, 0.2f, 1),
        v0, v1, v2, &h);
    std::printf("seg-tri hit=%d t=%.4f bary=(%.3f,%.3f,%.3f)\n",
                hit, h.t, h.bary.x, h.bary.y, h.bary.z);
    check(hit && std::fabs(h.t - 0.5f) < 1e-5f &&
          std::fabs(h.bary.x - 0.6f) < 1e-4f &&
          std::fabs(h.bary.y - 0.2f) < 1e-4f &&
          std::fabs(h.bary.z - 0.2f) < 1e-4f,
          "crossing segment has exact t and barycentrics");

    hit = geomlib_segment_triangle_intersect(
        geomlib_vec3_create(0.2f, 0.2f, -1), geomlib_vec3_create(0.2f, 0.2f, -0.1f),
        v0, v1, v2, &h);
    check(!hit, "segment stopping before the triangle misses");

    hit = geomlib_segment_triangle_intersect(
        geomlib_vec3_create(0.8f, 0.8f, -1), geomlib_vec3_create(0.8f, 0.8f, 1),
        v0, v1, v2, &h);
    check(!hit, "segment outside the triangle misses");

    hit = geomlib_segment_triangle_intersect(
        geomlib_vec3_create(-1, 0.2f, 0), geomlib_vec3_create(1, 0.2f, 0),
        v0, v1, v2, &h);
    check(!hit, "coplanar segment is not reported");
    if (failures != 0) {
        std::printf("%d segment-triangle test(s) failed\n", failures);
        return 1;
    }
    std::printf("All segment-triangle tests passed\n");
    return 0;
}
