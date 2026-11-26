/*
 * Anisotropic metrics: Euclidean limit, weighted distance, frame rotation,
 * determinant/anisotropy, curvature construction and mesh field.
 */
#include <geomlib/geomlib_metric.h>

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
    geomlib_vec3 o = geomlib_vec3_create(0, 0, 0);
    geomlib_vec3 q = geomlib_vec3_create(2, 3, 6);
    geomlib_metric3 I = geomlib_metric_identity();
    double d = geomlib_metric_distance(I, o, q);
    std::printf("metric identity d=%.6f expect=7\n", d);
    check(std::fabs(d - 7.0) < 1e-9, "identity metric is Euclidean");

    geomlib_metric3 D = geomlib_metric_diagonal(4.0, 1.0, 0.25);
    d = geomlib_metric_distance(D, o, q);
    std::printf("metric diagonal d=%.6f expect=%.6f\n", d, std::sqrt(34.0));
    check(std::fabs(d - std::sqrt(34.0)) < 1e-9, "diagonal metric weights axes");
    check(std::fabs(geomlib_metric_det(D) - 1.0) < 1e-9, "determinant is 1");
    check(std::fabs(geomlib_metric_anisotropy(D) - 4.0) < 1e-9,
          "anisotropy ratio is sqrt(16/1)");

    /* Rotated frame: rotating both the frame and the vector keeps distance. */
    geomlib_vec3 axes[3] = {
        geomlib_vec3_create(0, 1, 0),
        geomlib_vec3_create(-1, 0, 0),
        geomlib_vec3_create(0, 0, 1)
    };
    double eig[3] = {4.0, 1.0, 0.25};
    geomlib_metric3 R = geomlib_metric_from_frame(axes, eig);
    geomlib_vec3 v = geomlib_vec3_create(3, 2, 6);
    check(std::fabs(geomlib_metric_norm(R, v) -
                    geomlib_metric_norm(D, geomlib_vec3_create(2, -3, 6))) < 1e-9,
          "frame rotation preserves metric distance");

    geomlib_metric3 C;
    int ok = geomlib_metric_from_curvature(
        geomlib_vec3_create(0, 0, 1), 2.0f, 1.0f,
        geomlib_vec3_create(1, 0, 0), geomlib_vec3_create(0, 1, 0),
        0.1, 0.2, &C);
    std::printf("curvature metric aniso=%.6f\n", geomlib_metric_anisotropy(C));
    check(ok && geomlib_metric_anisotropy(C) > 1.0,
          "curvature metric is anisotropic when normal scale differs");

    geomlib_trimesh sphere;
    geomlib_create_sphere(&sphere, 1.0f, 16, 8);
    std::vector<geomlib_metric3> field(sphere.vertex_count);
    ok = geomlib_metric_field_from_curvature(&sphere, 0.2, 0.5, field.data());
    check(ok && geomlib_metric_det(field[0]) > 0.0,
          "mesh metric field is built for a sphere");
    geomlib_trimesh_free(&sphere);

    if (failures != 0) {
        std::printf("%d metric test(s) failed\n", failures);
        return 1;
    }
    std::printf("All metric tests passed\n");
    return 0;
}
