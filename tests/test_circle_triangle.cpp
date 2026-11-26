/*
 * Circle-triangle: interior centre, disjoint circle, containment and tangency.
 */
#include <geomlib/geomlib_circle_triangle.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec2 v0 = geomlib_vec2_create(0, 0);
    geomlib_vec2 v1 = geomlib_vec2_create(2, 0);
    geomlib_vec2 v2 = geomlib_vec2_create(0, 2);
    geomlib_circle_triangle_result r;

    int hit = geomlib_circle_triangle_intersect(
        geomlib_vec2_create(0.5f, 0.5f), 0.1f, v0, v1, v2, &r);
    std::printf("circle-tri interior hit=%d inside=%d sd=%.4f\n",
                hit, r.center_inside, r.signed_distance);
    check(hit && r.center_inside && r.signed_distance > 0.0f,
          "centre inside the triangle intersects even with positive edge gap");

    hit = geomlib_circle_triangle_intersect(
        geomlib_vec2_create(3, 3), 0.5f, v0, v1, v2, &r);
    check(!hit, "distant circle misses");

    hit = geomlib_circle_triangle_intersect(
        geomlib_vec2_create(0, 0), 2.0f,
        geomlib_vec2_create(1, 0), geomlib_vec2_create(0, 1),
        geomlib_vec2_create(1, 1), &r);
    check(hit && r.contains_triangle, "large circle contains the triangle");

    hit = geomlib_circle_triangle_intersect(
        geomlib_vec2_create(1, -1), 1.0f, v0, v1, v2, &r);
    std::printf("tangent sd=%.6f closest=(%.2f,%.2f)\n",
                r.signed_distance, r.closest.x, r.closest.y);
    check(hit && std::fabs(r.signed_distance) < 1e-5f &&
          std::fabs(r.closest.y) < 1e-5f,
          "circle tangent to an edge at (1,0)");
    if (failures != 0) {
        std::printf("%d circle-triangle test(s) failed\n", failures);
        return 1;
    }
    std::printf("All circle-triangle tests passed\n");
    return 0;
}
