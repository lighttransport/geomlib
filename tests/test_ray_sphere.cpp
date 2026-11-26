/*
 * Ray-sphere: entry/exit roots, misses, tangency and an interior origin.
 */
#include <geomlib/geomlib_ray_sphere.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_vec3 o = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 d = geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 c = geomlib_vec3_create(5, 0, 0);
    geomlib_ray_sphere_hit h;

    int hit = geomlib_ray_sphere_intersect(o, d, c, 1.0f, &h);
    std::printf("ray-sphere t0=%.4f t1=%.4f\n", h.t0, h.t1);
    check(hit && std::fabs(h.t0 - 4.0) < 1e-6 && std::fabs(h.t1 - 6.0) < 1e-6,
          "entry and exit roots are exact");
    check(std::fabs(h.p0.x - 4.0f) < 1e-5f &&
          std::fabs(h.n0.x + 1.0f) < 1e-5f &&
          std::fabs(h.n1.x - 1.0f) < 1e-5f,
          "points and outward normals are correct");

    hit = geomlib_ray_sphere_intersect(o, geomlib_vec3_create(0, 1, 0), c, 1.0f, &h);
    check(!hit, "ray misses sphere");

    /* Tangent from origin to sphere centred (0,5,0) radius 5 at t=5. */
    hit = geomlib_ray_sphere_intersect(o, d, geomlib_vec3_create(0, 5, 0), 5.0f, &h);
    check(hit && std::fabs(h.t0 - h.t1) < 1e-4,
          "tangent ray has coincident roots");

    /* Origin inside: near root negative, far root positive. */
    hit = geomlib_ray_sphere_intersect(c, d, c, 1.0f, &h);
    check(hit && h.t0 < 0.0 && h.t1 > 0.0 &&
          std::fabs(h.t0 + 1.0) < 1e-6 && std::fabs(h.t1 - 1.0) < 1e-6,
          "interior origin gives one negative and one positive root");

    if (failures != 0) {
        std::printf("%d ray-sphere test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ray-sphere tests passed\n");
    return 0;
}
