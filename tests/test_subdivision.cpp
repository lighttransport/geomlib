/* Subdivision: modified Butterfly stencils and Loop limit evaluation. */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static geomlib_trimesh make_plane(int n, float (*height)(int, int) = nullptr) {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, n * n, (n - 1) * (n - 1) * 2);
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++)
            geomlib_trimesh_add_vertex(&mesh,
                {(float)x, (float)y, height ? height(x, y) : 0.0f}, nullptr, nullptr);
    for (int y = 0; y + 1 < n; y++)
        for (int x = 0; x + 1 < n; x++) {
            uint32_t a = y * n + x, b = a + 1, c = a + n, d = c + 1;
            geomlib_trimesh_add_triangle(&mesh, a, b, d);
            geomlib_trimesh_add_triangle(&mesh, a, d, c);
        }
    return mesh;
}

static double mesh_area(const geomlib_trimesh *mesh) {
    double area = 0.0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(mesh, i, &a, &b, &c);
        area += 0.5 * geomlib_vec3_length(
            geomlib_vec3_cross(geomlib_vec3_sub(b, a), geomlib_vec3_sub(c, a)));
    }
    return area;
}

static int test_butterfly_plane() {
    geomlib_trimesh input = make_plane(7);
    std::vector<geomlib_vec3> control(input.vertices, input.vertices + input.vertex_count);
    double area = mesh_area(&input);
    geomlib_trimesh output = {};
    CHECK(geomlib_subdivide_butterfly(&input, &output) == GEOMLIB_MESH_SUCCESS,
          "Butterfly subdivision should succeed");
    CHECK(output.vertex_count > input.vertex_count && output.triangle_count == input.triangle_count * 4,
          "Butterfly splits every triangle into four");
    for (size_t i = 0; i < control.size(); i++)
        CHECK(std::memcmp(&control[i], &output.vertices[i], sizeof(geomlib_vec3)) == 0,
              "Butterfly is interpolating");
    for (size_t i = 0; i < output.vertex_count; i++)
        CHECK(std::fabs(output.vertices[i].z) < 1e-6f, "A planar grid stays planar");
    CHECK(std::fabs(mesh_area(&output) - area) < 1e-4 * area, "Butterfly preserves area");
    geomlib_trimesh_free(&output);
    geomlib_trimesh_free(&input);
    return 0;
}

static float curved_height(int x, int y) {
    (void)y;
    return 0.1f * (float)(x * x);
}

static int test_butterfly_stencil() {
    geomlib_trimesh input = make_plane(7, curved_height);
    geomlib_trimesh output = {};
    CHECK(geomlib_subdivide_butterfly(&input, &output) == GEOMLIB_MESH_SUCCESS,
          "Butterfly on a curved grid should succeed");
    /* At least one new edge vertex must differ from the linear midpoint, which
     * rules out plain midpoint subdivision. */
    bool nonlinear = false;
    for (size_t i = input.vertex_count; i < output.vertex_count && !nonlinear; i++) {
        float x = output.vertices[i].x, y = output.vertices[i].y;
        float fx = std::floor(x), fy = std::floor(y);
        int x0 = 0, y0 = 0, x1 = 0, y1 = 0;
        if (std::fabs(x - fx) < 1e-5f && std::fabs(y - fy - 0.5f) < 1e-5f) {
            x0 = (int)fx; y0 = (int)fy; x1 = (int)fx; y1 = (int)fy + 1;
        } else if (std::fabs(y - fy) < 1e-5f && std::fabs(x - fx - 0.5f) < 1e-5f) {
            x0 = (int)fx; y0 = (int)fy; x1 = (int)fx + 1; y1 = (int)fy;
        } else if (std::fabs(x - fx - 0.5f) < 1e-5f && std::fabs(y - fy - 0.5f) < 1e-5f) {
            x0 = (int)fx; y0 = (int)fy; x1 = (int)fx + 1; y1 = (int)fy + 1;
        } else {
            continue;
        }
        if (x0 < 0 || y0 < 0 || x1 >= 7 || y1 >= 7) continue;
        float midpoint = 0.5f * (curved_height(x0, y0) + curved_height(x1, y1));
        nonlinear |= std::fabs(output.vertices[i].z - midpoint) > 1e-3f;
    }
    CHECK(nonlinear, "Butterfly uses a non-linear edge stencil");
    geomlib_trimesh_free(&output);
    geomlib_trimesh_free(&input);
    return 0;
}

static geomlib_trimesh make_octahedron() {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, 6, 8);
    for (geomlib_vec3 p : {geomlib_vec3{1, 0, 0}, geomlib_vec3{0, 1, 0},
                           geomlib_vec3{-1, 0, 0}, geomlib_vec3{0, -1, 0},
                           geomlib_vec3{0, 0, 1}, geomlib_vec3{0, 0, -1}})
        geomlib_trimesh_add_vertex(&mesh, p, nullptr, nullptr);
    for (uint32_t i = 0; i < 4; i++) {
        geomlib_trimesh_add_triangle(&mesh, i, (i + 1) % 4, 4);
        geomlib_trimesh_add_triangle(&mesh, (i + 1) % 4, i, 5);
    }
    return mesh;
}

static int test_loop_limit() {
    geomlib_trimesh plane = make_plane(7);
    geomlib_vec3 interior = {3, 3, 0};
    geomlib_vec3 limit = geomlib_loop_limit_position(&plane, 3 * 7 + 3);
    CHECK(geomlib_vec3_length(geomlib_vec3_sub(limit, interior)) < 1e-4f,
          "The Loop limit of a symmetric planar vertex is the vertex");
    geomlib_vec3 normal = geomlib_loop_limit_normal(&plane, 3 * 7 + 3);
    CHECK(std::fabs(normal.z) > 0.99f, "The Loop limit normal follows the plane");
    geomlib_trimesh_free(&plane);

    geomlib_trimesh mesh = make_octahedron();
    geomlib_vec3 previous = mesh.vertices[0];
    double previous_distance = 1e30;
    for (int level = 0; level < 4; level++) {
        geomlib_vec3 limit0 = geomlib_loop_limit_position(&mesh, 0);
        double distance = geomlib_vec3_length(geomlib_vec3_sub(mesh.vertices[0], limit0));
        CHECK(distance < previous_distance + 1e-6,
              "Subdivision moves control vertices toward the Loop limit");
        previous_distance = distance;
        previous = mesh.vertices[0];
        (void)previous;
        geomlib_trimesh next = {};
        CHECK(geomlib_subdivide_loop(&mesh, &next) == GEOMLIB_MESH_SUCCESS,
              "Loop subdivision should succeed");
        geomlib_trimesh_free(&mesh);
        mesh = next;
    }
    geomlib_vec3 limit_final = geomlib_loop_limit_position(&mesh, 0);
    CHECK(geomlib_vec3_length(limit_final) > 0.0f, "Limit position stays finite");
    geomlib_trimesh_free(&mesh);
    return 0;
}

static int test_loop_boundary_and_projection() {
    geomlib_trimesh plane = make_plane(5);
    plane.vertices[0].z = 0.5f;
    geomlib_subdivision_params params =
        geomlib_subdivision_params_default(GEOMLIB_SUBDIV_LOOP, 1);
    params.preserve_boundaries = true;

    geomlib_trimesh frozen = {};
    CHECK(geomlib_subdivide_mesh(&plane, &params, &frozen) == GEOMLIB_MESH_SUCCESS,
          "Frozen-boundary Loop subdivision succeeds");
    CHECK(std::fabs(frozen.vertices[0].z - 0.5f) < 1e-6f,
          "preserve_boundaries keeps boundary vertices fixed");

    params.preserve_boundaries = false;
    geomlib_trimesh open = {};
    CHECK(geomlib_subdivide_mesh(&plane, &params, &open) == GEOMLIB_MESH_SUCCESS,
          "Boundary-rule Loop subdivision succeeds");
    CHECK(std::fabs(open.vertices[0].z - 0.375f) < 1e-6f,
          "Boundary vertices use 3/4 v + 1/8(n1+n2)");
    bool midpoint_found = false;
    for (size_t i = plane.vertex_count; i < open.vertex_count && !midpoint_found; i++)
        if (std::fabs(open.vertices[i].x - 0.5f) < 1e-6f &&
            std::fabs(open.vertices[i].y) < 1e-6f &&
            std::fabs(open.vertices[i].z - 0.25f) < 1e-6f)
            midpoint_found = true;
    CHECK(midpoint_found, "Boundary edges use their midpoint");
    geomlib_trimesh_free(&frozen);
    geomlib_trimesh_free(&open);
    geomlib_trimesh_free(&plane);

    geomlib_trimesh octahedron = make_octahedron();
    params = geomlib_subdivision_params_default(GEOMLIB_SUBDIV_LOOP, 2);
    params.project_to_limit = false;
    geomlib_trimesh subdivided = {};
    CHECK(geomlib_subdivide_mesh(&octahedron, &params, &subdivided) == GEOMLIB_MESH_SUCCESS,
          "Two-level Loop subdivision succeeds");
    params.project_to_limit = true;
    geomlib_trimesh projected = {};
    CHECK(geomlib_subdivide_mesh(&octahedron, &params, &projected) == GEOMLIB_MESH_SUCCESS,
          "Limit projection succeeds");
    CHECK(projected.vertex_count == subdivided.vertex_count, "Projection preserves topology");
    for (uint32_t i = 0; i < projected.vertex_count; i += 17) {
        geomlib_vec3 limit = geomlib_loop_limit_position(&subdivided, i);
        CHECK(geomlib_vec3_length(geomlib_vec3_sub(limit, projected.vertices[i])) < 1e-4f,
              "project_to_limit matches the Loop evaluation mask");
    }
    geomlib_trimesh_free(&subdivided);
    geomlib_trimesh_free(&projected);
    geomlib_trimesh_free(&octahedron);
    return 0;
}

static int test_catmull_clark() {
    /* create_box emits duplicated vertices per face; weld it first. */
    geomlib_trimesh raw_box = {};
    geomlib_create_box(&raw_box, 1.0f, 1.0f, 1.0f);
    geomlib_trimesh box = raw_box;
    geomlib_weld_result weld = geomlib_weld_vertices(&box, 1e-5f);
    geomlib_apply_weld_result(&box, &weld);
    geomlib_weld_result_free(&weld);
    geomlib_trimesh out = {};
    CHECK(geomlib_subdivide_catmull_clark(&box, &out) == GEOMLIB_MESH_SUCCESS,
          "Catmull-Clark subdivision should succeed");
    CHECK(out.triangle_count == box.triangle_count * 6,
          "Catmull-Clark splits each triangle into six triangles");
    CHECK(out.vertex_count > box.vertex_count,
          "Catmull-Clark adds face and edge points");
    /* The subdivided box stays closed and manifold. */
    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&out);
    if (!(analysis.is_closed && analysis.is_manifold))
        std::printf("  CC analysis: closed=%d manifold=%d edge=%d vertex=%d oriented=%d bad=%zu\n",
                    analysis.is_closed ? 1 : 0, analysis.is_manifold ? 1 : 0,
                    analysis.is_edge_manifold ? 1 : 0,
                    analysis.is_vertex_manifold ? 1 : 0,
                    analysis.is_oriented ? 1 : 0, analysis.nonmanifold_count);
    CHECK(analysis.is_closed && analysis.is_manifold,
          "Catmull-Clark preserves a closed manifold");
    geomlib_manifold_analysis_free(&analysis);
    geomlib_trimesh_free(&out);
    geomlib_trimesh_free(&box);
    return 0;
}

int main() {
    std::printf("=== Subdivision Test Suite ===\n\n");
    int failures = 0;
    failures += test_butterfly_plane();
    failures += test_butterfly_stencil();
    failures += test_loop_limit();
    failures += test_loop_boundary_and_projection();
    failures += test_catmull_clark();
    if (!failures) std::printf("=== ALL TESTS PASSED ===\n");
    else std::printf("=== %d TESTS FAILED ===\n", failures);
    return failures;
}
