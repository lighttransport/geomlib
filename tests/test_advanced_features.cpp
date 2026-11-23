/*
 * Comprehensive test for advanced geometry processing features
 */

#include <cstdio>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include "geomlib/geomlib.h"

#define TEST_ASSERT(cond, msg) \
    if (!(cond)) { \
        printf("FAIL: %s (line %d)\n", msg, __LINE__); \
        return 1; \
    }

static int test_hole_repair() {
    const geomlib_vec3 rim[] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {1, 1, 0}, {0, 2, 0}};
    for (int method = 0; method <= GEOMLIB_FILL_NATURAL; method++) {
        geomlib_trimesh mesh;
        geomlib_trimesh_init(&mesh, 0, 0);
        for (int i = 0; i < 5; i++) {
            geomlib_vec3 point = rim[i];
            if (method == GEOMLIB_FILL_NATURAL && i == 2) point.z = 0.35f;
            if (method == GEOMLIB_FILL_TANGENT_PLANE && i == 3) point.y = 2.4f;
            geomlib_trimesh_add_vertex(&mesh, point, nullptr, nullptr);
        }
        geomlib_trimesh_add_vertex(&mesh, {1, 0.7f, -2}, nullptr, nullptr);
        for (uint32_t i = 0; i < 5; i++) geomlib_trimesh_add_triangle(&mesh, i, 5, (i + 1) % 5);
        geomlib_vec3 original[6];
        std::memcpy(original, mesh.vertices, sizeof(original));
        geomlib_holes holes = geomlib_detect_holes(&mesh);
        TEST_ASSERT(holes.hole_count == 1 && holes.holes[0].vertex_count == 5,
                    "A five-vertex rim should form one closed hole");

        auto params = geomlib_fill_params_default(static_cast<geomlib_fill_method>(method));
        uint32_t invalid_loop[] = {0, 1, 5};
        geomlib_hole invalid = holes.holes[0];
        invalid.boundary_vertices = invalid_loop;
        invalid.vertex_count = 3;
        TEST_ASSERT(geomlib_fill_hole(&mesh, &invalid, &params) != GEOMLIB_MESH_SUCCESS &&
                    mesh.triangle_count == 5, "Interior edges cannot be used as a hole rim");

        TEST_ASSERT(geomlib_fill_hole(&mesh, &holes.holes[0], &params) == GEOMLIB_MESH_SUCCESS,
                    "Cap should be inserted");
        TEST_ASSERT(std::memcmp(original, mesh.vertices, sizeof(original)) == 0,
                    "Fairing must hold every original vertex fixed");
        geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&mesh);
        TEST_ASSERT(analysis.is_closed && analysis.is_manifold && analysis.is_oriented,
                    "Cap should close the mesh with consistent winding");
        if (method != GEOMLIB_FILL_TANGENT_PLANE) {
            TEST_ASSERT(mesh.vertex_count > 6, "Fairing methods should insert an interior");
        }
        geomlib_manifold_analysis_free(&analysis);
        geomlib_holes_free(&holes);
        geomlib_trimesh_free(&mesh);
    }

    const geomlib_vec2 clockwise[] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    uint32_t triangles[6];
    size_t count;
    TEST_ASSERT(geomlib_triangulate_polygon_2d(clockwise, 4, triangles, &count) == GEOMLIB_MESH_SUCCESS &&
                count == 2, "Public polygon triangulation should support clockwise winding");
    const geomlib_vec2 crossed[] = {{0, 0}, {1, 1}, {1, 0}, {0, 1}};
    TEST_ASSERT(geomlib_triangulate_polygon_2d(crossed, 4, triangles, &count) ==
                GEOMLIB_MESH_ERROR_DEGENERATE && count == 0, "Crossed polygons should fail cleanly");
    return 0;
}

static geomlib_sdf test_grid(float *values, uint32_t x, uint32_t y, uint32_t z) {
    geomlib_sdf sdf = {};
    sdf.values = values;
    sdf.voxel_size = 1.0f;
    sdf.resolution[0] = x;
    sdf.resolution[1] = y;
    sdf.resolution[2] = z;
    sdf.total_voxels = static_cast<size_t>(x) * y * z;
    return sdf;
}

static int test_marching_cubes() {
    float values[8];
    geomlib_sdf sdf = test_grid(values, 2, 2, 2);
    for (unsigned configuration = 0; configuration < 256; configuration++) {
        for (unsigned i = 0; i < 8; i++) {
            float magnitude = 1.0f + i * 0.19f;
            values[i] = (configuration & (1u << i)) ? -magnitude : magnitude;
        }
        geomlib_trimesh mesh;
        TEST_ASSERT(geomlib_sdf_to_mesh(&sdf, 0, &mesh) == GEOMLIB_MESH_SUCCESS,
                    "Every cube configuration should extract");
        TEST_ASSERT((mesh.triangle_count > 0) == (configuration != 0 && configuration != 255),
                    "Every mixed cube configuration should produce triangles");
        for (size_t i = 0; i < mesh.triangle_count * 3; i++)
            TEST_ASSERT(mesh.indices[i] < mesh.vertex_count, "Triangle vertices must exist");
        geomlib_trimesh_free(&mesh);
    }

    float ambiguous[12] = {1, -2, 1, 1, 1, 1, 1, 3, 1, 1, -4, 1};
    sdf = test_grid(ambiguous, 3, 2, 2);
    geomlib_trimesh mesh;
    TEST_ASSERT(geomlib_sdf_extract_surface(&sdf, &mesh) == GEOMLIB_MESH_SUCCESS,
                "Ambiguous neighbors should extract");
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(&mesh);
    size_t seams = 0;
    for (size_t i = 0; i < edges.edge_count; i++) {
        const auto &edge = edges.edges[i];
        if (mesh.vertices[edge.v0].x == 1.5f && mesh.vertices[edge.v1].x == 1.5f) {
            seams++;
            TEST_ASSERT(edge.face1 != UINT32_MAX, "Neighboring cubes must share ambiguous face segments");
        }
    }
    TEST_ASSERT(seams == 2 && geomlib_is_consistently_oriented(&mesh),
                "The ambiguous face must have two consistently wound shared segments");
    geomlib_edge_mesh_free(&edges);
    geomlib_trimesh_free(&mesh);

    float sphere[12 * 12 * 12];
    sdf = test_grid(sphere, 12, 12, 12);
    sdf.origin = {-6, -6, -6};
    for (uint32_t z = 0; z < 12; z++) {
        for (uint32_t y = 0; y < 12; y++) {
            for (uint32_t x = 0; x < 12; x++) {
                geomlib_vec3 p = {x - 5.5f, y - 5.5f, z - 5.5f};
                sphere[geomlib_sdf_get_index(&sdf, x, y, z)] = geomlib_vec3_length(p) - 4.2f;
            }
        }
    }
    TEST_ASSERT(geomlib_sdf_extract_surface(&sdf, &mesh) == GEOMLIB_MESH_SUCCESS,
                "Analytic sphere should extract");
    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&mesh);
    TEST_ASSERT(analysis.is_closed && analysis.is_manifold && analysis.is_oriented &&
                mesh.triangle_count > 0, "Sphere extraction must be watertight and consistently wound");
    for (size_t face = 0; face < mesh.triangle_count; face++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(&mesh, face, &a, &b, &c);
        geomlib_vec3 n = geomlib_vec3_cross(geomlib_vec3_sub(b, a), geomlib_vec3_sub(c, a));
        TEST_ASSERT(geomlib_vec3_dot(a, n) > 0.0f, "Sphere triangles must point outward");
    }
    geomlib_manifold_analysis_free(&analysis);
    geomlib_trimesh_free(&mesh);

    for (unsigned i = 0; i < 8; i++) values[i] = (i & 1) ? 0.75f : -0.25f;
    sdf = test_grid(values, 2, 2, 2);
    sdf.origin = {10, 20, 30};
    TEST_ASSERT(geomlib_sdf_extract_surface(&sdf, &mesh) == GEOMLIB_MESH_SUCCESS,
                "Offset plane should extract");
    for (size_t i = 0; i < mesh.vertex_count; i++)
        TEST_ASSERT(std::abs(mesh.vertices[i].x - 10.75f) < 1e-6f,
                    "Extraction must use the voxel-center coordinate convention");
    TEST_ASSERT(geomlib_sdf_query_distance(&sdf, {10.75f, 21.0f, 31.0f}) == 0,
                "Query and extraction must agree on world coordinates");
    geomlib_trimesh_free(&mesh);

    auto cropped = geomlib_sdf_extract_region(&sdf, {{11, 20, 30}, {12, 22, 32}});
    TEST_ASSERT(cropped.values && cropped.resolution[0] == 1 && cropped.total_voxels == 4 &&
                cropped.origin.x == 11.0f && cropped.values[0] == 0.75f,
                "Region extraction should copy samples in place");
    geomlib_sdf_free(&cropped);
    auto resized = geomlib_sdf_resize(&sdf, 4);
    TEST_ASSERT(resized.values && resized.resolution[0] == 4 && resized.voxel_size == 0.5f &&
                geomlib_sdf_query_distance(&resized, {10.75f, 21.0f, 31.0f}) == 0,
                "Resizing should interpolate at the new voxel centers");
    geomlib_sdf_free(&resized);

    for (unsigned i = 0; i < 8; i++) values[i] = (i & 1) ? 0.0f : -1.0f;
    TEST_ASSERT(geomlib_sdf_extract_surface(&sdf, &mesh) == GEOMLIB_MESH_SUCCESS &&
                mesh.triangle_count > 0, "Exact sample hits should form a surface");
    for (size_t i = 0; i < mesh.vertex_count; i++)
        TEST_ASSERT(mesh.vertices[i].x == 11.5f, "Exact iso hits should stay on the sample plane");
    geomlib_trimesh_free(&mesh);
    geomlib_trimesh tetra;
    geomlib_trimesh_init(&tetra, 4, 4);
    for (geomlib_vec3 p : {geomlib_vec3{0, 0, 0}, geomlib_vec3{1, 0, 0},
                         geomlib_vec3{0, 1, 0}, geomlib_vec3{0, 0, 1}}) {
        geomlib_trimesh_add_vertex(&tetra, p, nullptr, nullptr);
    }
    geomlib_trimesh_add_triangle(&tetra, 0, 2, 1);
    geomlib_trimesh_add_triangle(&tetra, 0, 1, 3);
    geomlib_trimesh_add_triangle(&tetra, 1, 2, 3);
    geomlib_trimesh_add_triangle(&tetra, 2, 0, 3);
    TEST_ASSERT(std::abs(geomlib_compute_winding_number(&tetra, {0.1f, 0.1f, 0.1f}) - 1) < 1e-6f &&
                std::abs(geomlib_compute_winding_number(&tetra, {1, 1, 1})) < 1e-6f,
                "Winding number should distinguish the inside of a closed surface");
    float distance;
    auto nearest = geomlib_find_closest_point_on_mesh(&tetra, nullptr, {-1, 0, 0}, &distance);
    TEST_ASSERT(nearest.x == 0 && nearest.y == 0 && nearest.z == 0 && distance == 1,
                "Closest point fallback should return the nearest point and distance");
    geomlib_sdf generated = geomlib_mesh_to_sdf_resolution(&tetra, 10);
    TEST_ASSERT(generated.values && geomlib_sdf_query_distance(&generated, {0.2f, 0.2f, 0.2f}) < 0 &&
                geomlib_sdf_query_distance(&generated, {0.5f, 0.5f, 0.5f}) > 0,
                "Mesh-generated field should have the correct inside and outside signs");
    geomlib_sdf_free(&generated);
    geomlib_trimesh_free(&tetra);
    return 0;
}

/* Test BVH-accelerated self-intersection detection */
static int test_self_intersection() {
    printf("Testing BVH-accelerated self-intersection detection...\n");
    
    /* Create a mesh with self-intersection */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 8, 4);
    
    /* First triangle */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 1, 0), NULL, NULL);
    
    /* Second triangle (intersecting) */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.3f, -0.1f), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.3f, 0.1f), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.7f, 0), NULL, NULL);
    
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 3, 4, 5);
    
    /* Detect intersections */
    geomlib_self_intersection intersections = geomlib_detect_self_intersections(&mesh);
    
    printf("  Found %zu intersection pairs\n", intersections.pair_count);
    TEST_ASSERT(intersections.has_intersections && intersections.pair_count == 1 &&
                intersections.pairs[0].face0 == 0 && intersections.pairs[0].face1 == 1,
                "Should report the intersecting face pair");
    
    geomlib_self_intersection_free(&intersections);
    auto strict = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_STRICT);
    TEST_ASSERT(!strict.is_valid && strict.error_count == 1, "Strict validation should report overlap");
    geomlib_validation_result_free(&strict);
    mesh.vertices[3] = {0.2f, 0.1f, 0};
    mesh.vertices[4] = {0.8f, 0.1f, 0};
    mesh.vertices[5] = {0.5f, 0.5f, 0};
    intersections = geomlib_detect_self_intersections(&mesh);
    TEST_ASSERT(intersections.pair_count == 1, "Detect coplanar overlap");
    geomlib_self_intersection_free(&intersections);
    mesh.indices[3] = 1;
    mesh.indices[4] = 0;
    mesh.vertices[5] = {0.5f, -0.5f, 0};
    intersections = geomlib_detect_self_intersections(&mesh);
    TEST_ASSERT(!intersections.has_intersections, "A normal shared edge is not self-intersection");
    geomlib_self_intersection_free(&intersections);
    mesh.vertices[5] = {0.5f, 0.5f, 0};
    intersections = geomlib_detect_self_intersections(&mesh);
    TEST_ASSERT(intersections.has_intersections, "Overlap beyond the shared edge should be detected");
    geomlib_self_intersection_free(&intersections);
    geomlib_trimesh_free(&mesh);
    
    printf("  PASS\n");
    return 0;
}

static int test_nonmanifold_splitting() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 7, 4);
    geomlib_vec3 normal = {0, 0, 1};
    geomlib_vec2 uv = {0.25f, 0.75f};
    for (geomlib_vec3 p : {geomlib_vec3{0, 0, 0}, geomlib_vec3{1, 0, 0},
                         geomlib_vec3{1, 1, 0}, geomlib_vec3{0, 1, 0},
                         geomlib_vec3{-1, 0, 0}, geomlib_vec3{-1, -1, 0}, geomlib_vec3{0, -1, 0}}) {
        geomlib_trimesh_add_vertex(&mesh, p, &normal, &uv);
    }
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 4, 5);
    geomlib_trimesh_add_triangle(&mesh, 0, 2, 3);
    geomlib_trimesh_add_triangle(&mesh, 0, 5, 6);

    auto before = geomlib_analyze_manifold(&mesh);
    TEST_ASSERT(before.is_edge_manifold && !before.is_vertex_manifold &&
                before.nonmanifold_count == 1 && before.nonmanifold_elements[0].index == 0,
                "A disconnected fan is a vertex-only manifold defect");
    geomlib_manifold_analysis_free(&before);
    auto loops = geomlib_extract_boundary_loops(&mesh);
    TEST_ASSERT(loops.loop_count == 0, "Boundary walks through a junction cannot be returned as loops");
    geomlib_boundary_loops_free(&loops);

    TEST_ASSERT(geomlib_split_nonmanifold_vertices(&mesh) == 1 && mesh.vertex_count == 8,
                "Split one disconnected fan");
    TEST_ASSERT(mesh.indices[0] == mesh.indices[6] && mesh.indices[3] == mesh.indices[9] &&
                mesh.indices[0] != mesh.indices[3], "Each fan must keep its own shared vertex");
    TEST_ASSERT(mesh.normals[7].z == normal.z && mesh.texcoords[7].x == uv.x &&
                mesh.texcoords[7].y == uv.y, "Duplicated vertex must retain its attributes");
    TEST_ASSERT(geomlib_split_nonmanifold_vertices(&mesh) == 0, "Splitting should be idempotent");
    geomlib_trimesh_free(&mesh);
    return 0;
}

static int test_feature_curves() {
    const uint32_t vertices[][2] = {{3, 4}, {2, 3}, {3, 5}, {1, 2}, {12, 10}, {10, 11}, {11, 12},
                                     {20, 21}, {21, 22}};
    constexpr size_t count = sizeof(vertices) / sizeof(vertices[0]);
    geomlib_feature_edge edges[count] = {};
    for (size_t i = 0; i < count; i++) {
        edges[i].v0 = vertices[i][0];
        edges[i].v1 = vertices[i][1];
    }
    uint32_t **curves;
    size_t *lengths, curve_count;
    TEST_ASSERT(geomlib_extract_feature_curves(edges, count, &curves, &lengths, &curve_count) ==
                GEOMLIB_MESH_SUCCESS && curve_count == 5, "Expect three branches, a cycle, and an open chain");
    unsigned used[count] = {};
    int closed = 0, length_three = 0;
    for (size_t i = 0; i < curve_count; i++) {
        closed += curves[i][0] == curves[i][lengths[i] - 1];
        length_three += lengths[i] == 3;
        for (size_t j = 1; j < lengths[i]; j++) {
            uint32_t a = curves[i][j - 1], b = curves[i][j];
            size_t k = 0;
            for (; k < count; k++)
                if ((a == vertices[k][0] && b == vertices[k][1]) ||
                    (a == vertices[k][1] && b == vertices[k][0])) break;
            TEST_ASSERT(k < count, "Each curve segment should be an input edge");
            used[k]++;
        }
        free(curves[i]);
    }
    for (unsigned n : used) TEST_ASSERT(n == 1, "Each input feature edge must appear once");
    TEST_ASSERT(closed == 1 && length_three == 2, "Curves should stop at junctions");
    free(curves);
    free(lengths);
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 2);
    for (geomlib_vec3 p : {geomlib_vec3{0, 0, 0}, geomlib_vec3{1, 0, 0},
                         geomlib_vec3{0, 1, 0}, geomlib_vec3{0, -1, 0.1f}}) {
        geomlib_trimesh_add_vertex(&mesh, p, nullptr, nullptr);
    }
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 1, 0, 3);
    auto segments = geomlib_segment_planar_regions(&mesh, 15, 0.05f);
    TEST_ASSERT(segments.region_count == 2, "Plane distance should split similarly oriented triangles");
    geomlib_regions_free(&segments);
    segments = geomlib_segment_planar_regions(&mesh, 15, 0.15f);
    TEST_ASSERT(segments.region_count == 1 && segments.regions[0].triangle_count == 2,
                "Both triangles should fit the wider plane tolerance");
    geomlib_regions_free(&segments);
    geomlib_trimesh_free(&mesh);
    return 0;
}

static bool has_issue(const geomlib_validation_result &result, geomlib_validation_issue_type type) {
    for (size_t i = 0; i < result.issue_count; i++) {
        if (result.issues[i].type == type) return true;
    }
    return false;
}

static int test_validation_and_statistics() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 0, 0);
    auto stats = geomlib_compute_mesh_stats(&mesh, 4);
    TEST_ASSERT(stats.edges.median_length == 0 && stats.angles.avg_angle_deg == 0,
                "Empty mesh statistics should be finite zeros");
    geomlib_mesh_stats_free(&stats);

    for (geomlib_vec3 p : {geomlib_vec3{0, 0, 0}, geomlib_vec3{1, 0, 0},
                         geomlib_vec3{0, 1, 0}, geomlib_vec3{0, -1, 0}, geomlib_vec3{0, 0, NAN}}) {
        geomlib_trimesh_add_vertex(&mesh, p, nullptr, nullptr);
    }
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 3);
    auto basic = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_BASIC);
    auto standard = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_STANDARD);
    TEST_ASSERT(basic.is_valid && !standard.is_valid &&
                has_issue(standard, GEOMLIB_ISSUE_NONFINITE_VERTEX) &&
                has_issue(standard, GEOMLIB_ISSUE_UNUSED_VERTEX),
                "Coordinate and unused-vertex checks start at Standard");
    geomlib_validation_result_free(&basic);
    geomlib_validation_result_free(&standard);

    mesh.vertex_count--;
    standard = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_STANDARD);
    TEST_ASSERT(!standard.is_valid && has_issue(standard, GEOMLIB_ISSUE_INCONSISTENT_ORIENTATION),
                "Standard must detect inconsistent winding");
    geomlib_validation_result_free(&standard);
    TEST_ASSERT(geomlib_orient_consistently(&mesh) == GEOMLIB_MESH_SUCCESS &&
                geomlib_is_consistently_oriented(&mesh), "Orientation repair should correct the winding");

    stats = geomlib_compute_mesh_stats(&mesh, 4);
    TEST_ASSERT(stats.areas.bin_width == 0 && stats.areas.histogram[0] == 2 &&
                stats.edges.median_length == 1, "Uniform values should go in the first histogram bin");
    geomlib_mesh_stats_free(&stats);

    geomlib_trimesh_add_triangle(&mesh, 0, 1, 3);
    standard = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_STANDARD);
    TEST_ASSERT(has_issue(standard, GEOMLIB_ISSUE_NON_MANIFOLD_EDGE) &&
                has_issue(standard, GEOMLIB_ISSUE_NON_MANIFOLD_VERTEX), "Report nonmanifold edge and vertices");
    geomlib_validation_result_free(&standard);
    mesh.triangle_count = 2;
    mesh.indices[5] = UINT32_MAX;
    basic = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_BASIC);
    TEST_ASSERT(!basic.is_valid && has_issue(basic, GEOMLIB_ISSUE_INVALID_INDEX),
                "Basic validation must report out-of-range indices");
    geomlib_validation_result_free(&basic);
    geomlib_trimesh_free(&mesh);
    return 0;
}

/* Test surface reconstruction methods */
static int test_surface_reconstruction() {
    printf("Testing surface reconstruction (ball pivoting & alpha shape)...\n");
    
    /* Create a simple point cloud (cube corners) */
    geomlib_pointcloud pc;
    pc.point_count = 8;
    pc.points = (geomlib_vec3*)malloc(8 * sizeof(geomlib_vec3));
    pc.normals = NULL;
    
    /* Cube corners */
    pc.points[0] = geomlib_vec3_create(0, 0, 0);
    pc.points[1] = geomlib_vec3_create(1, 0, 0);
    pc.points[2] = geomlib_vec3_create(0, 1, 0);
    pc.points[3] = geomlib_vec3_create(1, 1, 0);
    pc.points[4] = geomlib_vec3_create(0, 0, 1);
    pc.points[5] = geomlib_vec3_create(1, 0, 1);
    pc.points[6] = geomlib_vec3_create(0, 1, 1);
    pc.points[7] = geomlib_vec3_create(1, 1, 1);
    
    /* Test ball pivoting */
    geomlib_surface_params params = geomlib_surface_params_default(GEOMLIB_SURFACE_BALL_PIVOTING);
    params.ball_radius = 1.5f;  /* Large enough to connect corners */
    
    geomlib_trimesh mesh;
    int result = geomlib_reconstruct_surface(&pc, &params, &mesh);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Ball pivoting should succeed");
    
    printf("  Ball pivoting: %zu triangles\n", mesh.triangle_count);
    geomlib_trimesh_free(&mesh);
    
    /* Test alpha shape */
    params = geomlib_surface_params_default(GEOMLIB_SURFACE_ALPHA_SHAPE);
    params.alpha = 1.5f;
    
    result = geomlib_reconstruct_surface(&pc, &params, &mesh);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Alpha shape should succeed");
    
    printf("  Alpha shape: %zu triangles\n", mesh.triangle_count);
    geomlib_trimesh_free(&mesh);
    
    free(pc.points);
    
    printf("  PASS\n");
    return 0;
}

/* Test SDF operations */
static int test_sdf_operations() {
    printf("Testing SDF operations (smoothing, blending, I/O)...\n");

    /* Create a simple mesh to convert to SDF */
    geomlib_trimesh mesh;
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);

    /* Create SDF from mesh */
    geomlib_sdf sdf1 = geomlib_mesh_to_sdf_resolution(&mesh, 16);
    TEST_ASSERT(sdf1.values != NULL, "SDF creation should succeed");

    geomlib_sdf sdf2 = geomlib_mesh_to_sdf_resolution(&mesh, 16);

    /* Test smoothing */
    int result = geomlib_sdf_smooth(&sdf1, 1.0f);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "SDF smoothing should succeed");

    /* Test blending */
    geomlib_sdf blended = geomlib_sdf_blend(&sdf1, &sdf2, 0.5f);
    TEST_ASSERT(blended.values != NULL, "SDF blending should succeed");

    /* Test file I/O */
    result = geomlib_sdf_save(&sdf1, "/tmp/test_sdf.bin");
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "SDF save should succeed");

    geomlib_sdf loaded;
    result = geomlib_sdf_load(&loaded, "/tmp/test_sdf.bin");
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "SDF load should succeed");
    TEST_ASSERT(loaded.resolution[0] == sdf1.resolution[0], "Resolution should match");

    geomlib_trimesh_free(&mesh);
    geomlib_sdf_free(&sdf1);
    geomlib_sdf_free(&sdf2);
    geomlib_sdf_free(&blended);
    geomlib_sdf_free(&loaded);

    printf("  PASS\n");
    return 0;
}

static int test_pointcloud_preprocessing() {
    geomlib_pointcloud cloud;
    TEST_ASSERT(geomlib_pointcloud_init(&cloud, 0) == GEOMLIB_MESH_SUCCESS,
                "Empty point cloud should be appendable");
    for (int y = -2; y <= 2; y++) {
        for (int x = -2; x <= 2; x++) {
            float z = float(2 * x - y);
            geomlib_pointcloud_add_point(&cloud, {float(x), float(y), z}, nullptr, 1);
        }
    }
    TEST_ASSERT(geomlib_pointcloud_estimate_normals(&cloud, 6) == GEOMLIB_MESH_SUCCESS,
                "Estimate PCA normals on a tilted plane");
    geomlib_vec3 expected = geomlib_vec3_normalize({-2, 1, 1});
    for (size_t i = 0; i < cloud.point_count; i++) {
        TEST_ASSERT(std::abs(geomlib_vec3_dot(expected, cloud.normals[i])) > 0.9999f &&
                    geomlib_vec3_dot(cloud.normals[0], cloud.normals[i]) > 0.9999f,
                    "Plane normals have the correct direction and one consistent sign");
    }
    geomlib_pointcloud_free(&cloud);

    geomlib_pointcloud_init(&cloud, 4);
    geomlib_vec3 normal = {0, 1, 0};
    geomlib_pointcloud_add_point(&cloud, {-1, 0, 0}, &normal, 0.25f);
    geomlib_pointcloud_add_point(&cloud, {-0.5f, 0, 0}, &normal, 0.75f);
    geomlib_pointcloud_add_point(&cloud, {0.25f, 0, 0}, &normal, 0.1f);
    geomlib_pointcloud_add_point(&cloud, {0.75f, 0, 0}, &normal, 0.3f);
    geomlib_pointcloud sampled = {};
    TEST_ASSERT(geomlib_pointcloud_downsample(&cloud, 1.0f, &sampled) == GEOMLIB_MESH_SUCCESS &&
                sampled.point_count == 2, "Voxel downsampling should group points by relative cell");
    TEST_ASSERT(sampled.points[0].x == -0.75f && sampled.points[1].x == 0.5f &&
                sampled.normals[0].y == 1 && std::abs(sampled.confidences[1] - 0.2f) < 1e-6f,
                "Voxel averages should preserve positions, normals, and confidence");
    TEST_ASSERT(geomlib_pointcloud_estimate_normals(&cloud, 3) == GEOMLIB_MESH_SUCCESS &&
                geomlib_vec3_length(cloud.normals[0]) == 0, "Collinear neighborhoods have no normal");
    geomlib_pointcloud_free(&sampled);
    geomlib_pointcloud_free(&cloud);
    return 0;
}

static geomlib_trimesh make_hole_fixture(float raised_z) {
    const geomlib_vec3 rim[] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {1, 1, 0}, {0, 2, 0}};
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 0, 0);
    for (int i = 0; i < 5; i++) {
        geomlib_vec3 p = rim[i];
        if (i == 2) p.z = raised_z;
        geomlib_trimesh_add_vertex(&mesh, p, nullptr, nullptr);
    }
    geomlib_trimesh_add_vertex(&mesh, {1, 0.7f, -2}, nullptr, nullptr);
    for (uint32_t i = 0; i < 5; i++)
        geomlib_trimesh_add_triangle(&mesh, i, 5, (i + 1) % 5);
    return mesh;
}

static int test_hole_planarity() {
    geomlib_trimesh planar = make_hole_fixture(0.0f);
    geomlib_holes holes = geomlib_detect_holes(&planar);
    TEST_ASSERT(holes.hole_count == 1, "Planar fixture has one hole");
    const geomlib_hole &flat = holes.holes[0];
    TEST_ASSERT(flat.is_planar && flat.planarity > 0.99f,
                "A flat rim has high planarity");
    geomlib_vec3 origin, normal;
    geomlib_hole_best_fit_plane(&planar, &flat, &origin, &normal);
    TEST_ASSERT(std::fabs(normal.z) > 0.99f, "Best-fit plane matches the rim plane");
    geomlib_holes_free(&holes);
    geomlib_trimesh_free(&planar);

    geomlib_trimesh bent = make_hole_fixture(1.0f);
    holes = geomlib_detect_holes(&bent);
    TEST_ASSERT(holes.hole_count == 1, "Bent fixture has one hole");
    const geomlib_hole &curved = holes.holes[0];
    TEST_ASSERT(!curved.is_planar && curved.planarity < 0.9f,
                "A raised rim lowers the planarity measure");
    geomlib_holes_free(&holes);
    geomlib_trimesh_free(&bent);
    return 0;
}

int main() {
    printf("=== Advanced Features Test Suite ===\n\n");

    int failures = 0;

    failures += test_hole_repair();
    failures += test_hole_planarity();
    failures += test_marching_cubes();
    failures += test_self_intersection();
    failures += test_nonmanifold_splitting();
    failures += test_feature_curves();
    failures += test_validation_and_statistics();
    failures += test_surface_reconstruction();
    failures += test_pointcloud_preprocessing();
    failures += test_sdf_operations();
    
    printf("\n");
    if (failures == 0) {
        printf("=== ALL TESTS PASSED ===\n");
    } else {
        printf("=== %d TESTS FAILED ===\n", failures);
    }
    
    return failures;
}
