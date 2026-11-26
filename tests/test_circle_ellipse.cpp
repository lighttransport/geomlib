/*
 * Circle-ellipse: miss, penetration and centre containment.
 */
#include <geomlib/geomlib_circle_ellipse.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec2 o = geomlib_vec2_create(0, 0);
    geomlib_circle_ellipse_result r;

    int hit = geomlib_circle_ellipse_overlap(
        geomlib_vec2_create(3, 0), 0.5f, o, 2.0f, 1.0f, 0.0f, &r);
    std::printf("circle-ellipse miss d=%.4f closest=(%.2f,%.2f)\n",
                r.distance, r.closest.x, r.closest.y);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-3f &&
          std::fabs(r.closest.x - 2.0f) < 1e-3f,
          "circle one unit from the ellipse misses with radius 0.5");

    hit = geomlib_circle_ellipse_overlap(
        geomlib_vec2_create(3, 0), 1.5f, o, 2.0f, 1.0f, 0.0f, &r);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-3f,
          "larger circle penetrates by 0.5");

    hit = geomlib_circle_ellipse_overlap(o, 0.1f, o, 2.0f, 1.0f, 0.0f, &r);
    check(hit && r.center_inside && r.distance > 0.0f,
          "centre containment overlaps despite a positive boundary gap");
    if (failures != 0) {
        std::printf("%d circle-ellipse test(s) failed\n", failures);
        return 1;
    }
    std::printf("All circle-ellipse tests passed\n");
    return 0;
}
