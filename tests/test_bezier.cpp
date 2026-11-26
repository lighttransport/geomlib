/*
 * Bezier: endpoints are interpolated, linear controls reproduce the line,
 * subdivision endpoints join at the evaluated point, derivative is exact.
 */
#include <geomlib/geomlib_bezier.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    std::vector<geomlib_vec3> cp = {
        geomlib_vec3_create(0, 0, 0), geomlib_vec3_create(1, 0, 0),
        geomlib_vec3_create(2, 0, 0), geomlib_vec3_create(3, 0, 0)
    };
    geomlib_vec3 p0 = geomlib_bezier_evaluate(cp.data(), cp.size(), 0.0);
    geomlib_vec3 p1 = geomlib_bezier_evaluate(cp.data(), cp.size(), 1.0);
    geomlib_vec3 pm = geomlib_bezier_evaluate(cp.data(), cp.size(), 0.5);
    std::printf("bezier p0=(%.2f) p1=(%.2f) pm=(%.2f)\n", p0.x, p1.x, pm.x);
    check(std::fabs(p0.x) < 1e-6 && std::fabs(p1.x - 3.0f) < 1e-6,
          "endpoints interpolated");
    check(std::fabs(pm.x - 1.5f) < 1e-5, "linear controls give linear curve");

    geomlib_vec3 d0 = geomlib_bezier_derivative(cp.data(), cp.size(), 0.0);
    check(std::fabs(d0.x - 3.0f) < 1e-5, "derivative at start is 3(P1-P0)");

    std::vector<geomlib_vec3> left(cp.size()), right(cp.size());
    geomlib_bezier_subdivide(cp.data(), cp.size(), 0.5,
                             left.data(), right.data());
    check(geomlib_vec3_length(geomlib_vec3_sub(left[0], p0)) < 1e-5f &&
          geomlib_vec3_length(geomlib_vec3_sub(right.back(), p1)) < 1e-5f &&
          geomlib_vec3_length(geomlib_vec3_sub(left.back(), right[0])) < 1e-5f &&
          geomlib_vec3_length(geomlib_vec3_sub(left.back(), pm)) < 1e-5f,
          "subdivision halves join at the evaluated point");

    /* Re-evaluating the left half at t=1 must match the original at t=0.5. */
    for (double t = 0.0; t <= 1.0001; t += 0.25) {
        geomlib_vec3 a = geomlib_bezier_evaluate(left.data(), left.size(), t);
        geomlib_vec3 b = geomlib_bezier_evaluate(cp.data(), cp.size(), 0.5 * t);
        if (geomlib_vec3_length(geomlib_vec3_sub(a, b)) > 1e-4f) {
            check(false, "left half reproduces the original curve");
            break;
        }
    }

    if (failures != 0) {
        std::printf("%d Bezier test(s) failed\n", failures);
        return 1;
    }
    std::printf("All Bezier tests passed\n");
    return 0;
}
