/*
 * Circle-AABB: face projection, penetration and interior centre.
 */
#include <geomlib/geomlib_circle_aabb.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec2 lo = geomlib_vec2_create(0, 0);
    geomlib_vec2 hi = geomlib_vec2_create(1, 1);
    geomlib_circle_aabb_result r;

    int hit = geomlib_circle_aabb_overlap(
        geomlib_vec2_create(2, 0.5f), 0.5f, lo, hi, &r);
    std::printf("circle-aabb miss d=%.4f closest=(%.2f,%.2f)\n",
                r.distance, r.closest.x, r.closest.y);
    check(!hit && std::fabs(r.distance - 0.5f) < 1e-5f &&
          std::fabs(r.closest.x - 1.0f) < 1e-5f,
          "circle outside the rectangle misses");

    hit = geomlib_circle_aabb_overlap(
        geomlib_vec2_create(2, 0.5f), 1.5f, lo, hi, &r);
    check(hit && std::fabs(r.distance + 0.5f) < 1e-5f,
          "circle reaching the rectangle penetrates");

    hit = geomlib_circle_aabb_overlap(
        geomlib_vec2_create(0.5f, 0.5f), 0.2f, lo, hi, &r);
    check(hit && std::fabs(r.distance + 0.2f) < 1e-5f,
          "interior centre penetrates by the radius");
    if (failures != 0) {
        std::printf("%d circle-AABB test(s) failed\n", failures);
        return 1;
    }
    std::printf("All circle-AABB tests passed\n");
    return 0;
}
