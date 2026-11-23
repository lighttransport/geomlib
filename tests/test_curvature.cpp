/* Curvature: principal directions are tangent, orthonormal, and align with a cylinder axis. */
#include <cmath>
#include <cstdio>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static int check_tangent_frame(const geomlib_trimesh *mesh,
                                const geomlib_mesh_curvature *curvature, size_t stride) {
    for (size_t v = 0; v < mesh->vertex_count; v += stride) {
        geomlib_vec3 normal = mesh->normals[v];
        geomlib_vec3 dmax = curvature->curvatures[v].principal_dir_max;
        geomlib_vec3 dmin = curvature->curvatures[v].principal_dir_min;
        if (geomlib_vec3_length(normal) < 0.5f) continue;
        if (geomlib_vec3_length(dmax) < 0.5f || geomlib_vec3_length(dmin) < 0.5f) {
            std::printf("  zero direction at vertex %zu (normal %.3f)\n", v,
                        geomlib_vec3_length(normal));
            continue;
        }
        CHECK(std::fabs(geomlib_vec3_length(dmax) - 1.0f) < 1e-4f &&
              std::fabs(geomlib_vec3_length(dmin) - 1.0f) < 1e-4f,
              "Principal directions are unit length");
        CHECK(std::fabs(geomlib_vec3_dot(dmax, normal)) < 1e-3f &&
              std::fabs(geomlib_vec3_dot(dmin, normal)) < 1e-3f,
              "Principal directions lie in the tangent plane");
        CHECK(std::fabs(geomlib_vec3_dot(dmax, dmin)) < 1e-3f,
              "Principal directions are orthogonal");
    }
    return 0;
}

int main() {
    std::printf("=== Curvature Test Suite ===\n\n");

    geomlib_trimesh sphere = {};
    geomlib_create_sphere(&sphere, 1.0f, 16, 8);
    geomlib_trimesh_compute_normals(&sphere);
    geomlib_mesh_curvature sphere_curvature = geomlib_curvature_create(sphere.vertex_count);
    CHECK(geomlib_compute_curvature(&sphere, &sphere_curvature) == GEOMLIB_MESH_SUCCESS,
          "Sphere curvature succeeds");
    CHECK(check_tangent_frame(&sphere, &sphere_curvature, 7) == 0,
          "Sphere principal frame is valid");
    geomlib_curvature_free(&sphere_curvature);
    geomlib_trimesh_free(&sphere);

    geomlib_trimesh cylinder = {};
    geomlib_create_cylinder(&cylinder, 1.0f, 2.0f, 16);
    geomlib_trimesh_compute_normals(&cylinder);
    geomlib_mesh_curvature cylinder_curvature = geomlib_curvature_create(cylinder.vertex_count);
    CHECK(geomlib_compute_curvature(&cylinder, &cylinder_curvature) == GEOMLIB_MESH_SUCCESS,
          "Cylinder curvature succeeds");
    CHECK(check_tangent_frame(&cylinder, &cylinder_curvature, 11) == 0,
          "Cylinder principal frame is valid");

    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(&cylinder);
    float extent[3] = {bounds.max.x - bounds.min.x, bounds.max.y - bounds.min.y,
                       bounds.max.z - bounds.min.z};
    int axis = 0;
    for (int i = 1; i < 3; i++) if (extent[i] > extent[axis]) axis = i;
    geomlib_vec3 axis_dir = {axis == 0 ? 1.0f : 0.0f, axis == 1 ? 1.0f : 0.0f,
                             axis == 2 ? 1.0f : 0.0f};
    int checked = 0;
    for (size_t v = 0; v < cylinder.vertex_count; v++) {
        geomlib_vec3 p = cylinder.vertices[v];
        float height = (&p.x)[axis];
        float radial = 0.0f;
        for (int i = 0; i < 3; i++)
            if (i != axis) radial += (&p.x)[i] * (&p.x)[i];
        radial = std::sqrt(radial);
        if (radial < 0.5f || std::fabs(height) > 0.5f) continue;
        geomlib_vec3 dmax = cylinder_curvature.curvatures[v].principal_dir_max;
        geomlib_vec3 dmin = cylinder_curvature.curvatures[v].principal_dir_min;
        float align = std::fmax(std::fabs(geomlib_vec3_dot(dmax, axis_dir)),
                                std::fabs(geomlib_vec3_dot(dmin, axis_dir)));
        CHECK(align > 0.8f, "One principal direction follows the cylinder axis");
        checked++;
    }
    CHECK(checked > 8, "Cylinder side vertices were checked");
    std::printf("  checked %d cylinder side vertices\n", checked);

    geomlib_curvature_free(&cylinder_curvature);
    geomlib_trimesh_free(&cylinder);
    std::printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
