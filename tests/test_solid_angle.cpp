/*
 * Solid angle: an octant triangle subtends pi/2, a box gives 4pi inside and
 * 0 outside, and winding number matches.
 */
#include <geomlib/geomlib_solid_angle.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    const double PI = 3.14159265358979323846;
    geomlib_vec3 o = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 ax = geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 ay = geomlib_vec3_create(0, 1, 0);
    geomlib_vec3 az = geomlib_vec3_create(0, 0, 1);

    /* Triangle spanning one octant seen from the origin: 4pi/8 = pi/2. */
    double omega = geomlib_solid_angle_triangle(o, ax, ay, az);
    std::printf("octant solid angle=%.6f expect=%.6f\n", omega, PI / 2.0);
    check(std::fabs(omega - PI / 2.0) < 1e-8, "octant triangle subtends pi/2");

    /* Reversing winding flips the sign. */
    double rev = geomlib_solid_angle_triangle(o, az, ay, ax);
    check(std::fabs(rev + omega) < 1e-8, "reversing winding flips sign");

    /* A closed box: 4pi inside, 0 outside. */
    geomlib_trimesh box;
    geomlib_create_box(&box, 2.0f, 2.0f, 2.0f);
    double inside = geomlib_solid_angle_mesh(&box, o);
    double outside = geomlib_solid_angle_mesh(&box, geomlib_vec3_create(10, 10, 10));
    std::printf("box inside=%.6f (4pi=%.6f) outside=%.6f\n",
                inside, 4.0 * PI, outside);
    check(std::fabs(inside - 4.0 * PI) < 1e-3, "interior point subtends 4pi");
    check(std::fabs(outside) < 1e-3, "exterior point subtends 0");

    double wn = geomlib_winding_number_mesh(&box, o);
    check(std::fabs(wn - 1.0) < 1e-4, "winding number is 1 inside");
    wn = geomlib_winding_number_mesh(&box, geomlib_vec3_create(10, 10, 10));
    check(std::fabs(wn) < 1e-4, "winding number is 0 outside");
    geomlib_trimesh_free(&box);

    if (failures != 0) {
        std::printf("%d solid angle test(s) failed\n", failures);
        return 1;
    }
    std::printf("All solid angle tests passed\n");
    return 0;
}
