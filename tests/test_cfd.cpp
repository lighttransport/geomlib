/* Analytic field operators and tetrahedral finite-volume geometry. */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <fstream>
#include <string>
#include <vector>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static bool near(float a, float b) {
    return std::abs(a - b) < 1e-5f;
}

static bool near(geomlib_vec3 a, geomlib_vec3 b) {
    return near(a.x, b.x) && near(a.y, b.y) && near(a.z, b.z);
}

static int test_surface_operators() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 10, 8);
    for (int y = -1; y <= 1; y++)
        for (int x = -1; x <= 1; x++)
            geomlib_trimesh_add_vertex(&mesh, {float(x), float(y), 0}, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, {4, 4, 0}, nullptr, nullptr);
    for (uint32_t y = 0; y < 2; y++) {
        for (uint32_t x = 0; x < 2; x++) {
            uint32_t a = y * 3 + x;
            geomlib_trimesh_add_triangle(&mesh, a, a + 1, a + 4);
            geomlib_trimesh_add_triangle(&mesh, a, a + 4, a + 3);
        }
    }
    auto scalar = geomlib_scalar_field_create(mesh.vertex_count);
    auto scalar_out = geomlib_scalar_field_create(mesh.vertex_count);
    auto vector = geomlib_vector_field_create(mesh.vertex_count);
    auto vector_out = geomlib_vector_field_create(mesh.vertex_count);
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        auto p = mesh.vertices[i];
        scalar.values[i] = 7 + 2 * p.x - 3 * p.y;
        vector.values[i] = {2 * p.x - p.y, p.x + 3 * p.y, 5};
    }
    for (auto method : {GEOMLIB_GRADIENT_GREEN_GAUSS, GEOMLIB_GRADIENT_GREEN_GAUSS_NODE,
                       GEOMLIB_GRADIENT_LEAST_SQUARES}) {
        geomlib_vector_field gradient = {};
        CHECK(geomlib_compute_gradient_scalar(&mesh, &scalar, method, &gradient) ==
              GEOMLIB_MESH_SUCCESS, "Surface gradient should allocate its result");
        for (size_t i = 0; i < 9; i++)
            CHECK(near(gradient.values[i], {2, -3, 0}), "The affine tangential gradient is exact");
        CHECK(near(gradient.values[9], {}), "Isolated vertex has zero gradient");
        geomlib_vector_field_free(&gradient);
    }

    CHECK(geomlib_compute_divergence(&mesh, &vector, &scalar_out) == GEOMLIB_MESH_SUCCESS,
          "Divergence should succeed");
    CHECK(geomlib_compute_curl(&mesh, &vector, &vector_out) == GEOMLIB_MESH_SUCCESS,
          "Curl should succeed");
    for (size_t i = 0; i < 9; i++) {
        CHECK(near(scalar_out.values[i], 5), "Divergence is the tangential trace");
        CHECK(near(vector_out.values[i], {0, 0, 2}), "Curl has the analytic normal component");
    }
    CHECK(geomlib_reverse_orientation(&mesh) == GEOMLIB_MESH_SUCCESS, "Reverse fixture winding");
    CHECK(geomlib_compute_curl(&mesh, &vector, &vector) == GEOMLIB_MESH_SUCCESS &&
          near(vector.values[4], {0, 0, 2}), "Curl supports aliases and is invariant under winding reversal");

    CHECK(geomlib_compute_laplacian(&mesh, &scalar, &scalar_out) == GEOMLIB_MESH_SUCCESS &&
          near(scalar_out.values[4], 0), "An affine scalar field has zero interior Laplacian");
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        auto p = mesh.vertices[i];
        scalar.values[i] = p.x * p.x + p.y * p.y;
        vector.values[i] = {1, 2, 3};
    }
    CHECK(geomlib_compute_laplacian(&mesh, &scalar, &scalar) == GEOMLIB_MESH_SUCCESS &&
          near(scalar.values[4], 4) && near(scalar.values[9], 0),
          "Quadratic Laplace-Beltrami value is 4 at the interior vertex");
    CHECK(geomlib_compute_divergence(&mesh, &vector, &scalar_out) == GEOMLIB_MESH_SUCCESS &&
          geomlib_compute_curl(&mesh, &vector, &vector_out) == GEOMLIB_MESH_SUCCESS,
          "Constant vector operators should succeed");
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        CHECK(near(scalar_out.values[i], 0) && near(vector_out.values[i], {}),
              "Constant vector field has zero divergence and curl");
        scalar.values[i] = 12;
    }
    CHECK(geomlib_compute_laplacian(&mesh, &scalar, &scalar_out) == GEOMLIB_MESH_SUCCESS,
          "Constant scalar operator should succeed");
    for (size_t i = 0; i < mesh.vertex_count; i++)
        CHECK(near(scalar_out.values[i], 0), "Constant scalar field has zero Laplacian");

    mesh.indices[0] = UINT32_MAX;
    scalar_out.values[0] = 55;
    CHECK(geomlib_compute_divergence(&mesh, &vector, &scalar_out) == GEOMLIB_MESH_ERROR_INVALID &&
          scalar_out.values[0] == 55, "Bad index must leave output unchanged");
    geomlib_scalar_field_free(&scalar);
    geomlib_scalar_field_free(&scalar_out);
    geomlib_vector_field_free(&vector);
    geomlib_vector_field_free(&vector_out);
    geomlib_trimesh_free(&mesh);

    const geomlib_vec3 polygon[] = {{0, 0, 0}, {0, 1, 0}, {2, 1, 0}, {2, 0, 0}};
    const float values[] = {7, 4, 8, 11};
    CHECK(near(geomlib_gradient_green_gauss(polygon, values, 4), {2, -3, 0}),
          "Green-Gauss polygon flux recovers an affine field");
    return 0;
}

static int test_interpolation() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 6, 2);
    for (auto p : {geomlib_vec3{0, 0, 0}, geomlib_vec3{1, 0, 0}, geomlib_vec3{0, 1, 0},
                   geomlib_vec3{2, 0, 0}, geomlib_vec3{3, 0, 0}, geomlib_vec3{2, 1, 0}})
        geomlib_trimesh_add_vertex(&mesh, p, nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 3, 4, 5);
    auto scalar = geomlib_scalar_field_create(mesh.vertex_count);
    auto vector = geomlib_vector_field_create(mesh.vertex_count);
    for (size_t i = 3; i < mesh.vertex_count; i++) {
        auto p = mesh.vertices[i];
        scalar.values[i] = 1 + 2 * p.x + 4 * p.y;
        vector.values[i] = p;
    }
    CHECK(near(geomlib_interpolate_scalar(&mesh, &scalar, {2.25f, 0.25f, 0.3f},
                                         GEOMLIB_INTERP_LINEAR), 6.5f),
          "Scalar interpolation uses the closest face and its projected point");
    CHECK(near(geomlib_interpolate_vector(&mesh, &vector, {3, 1, 0},
                                         GEOMLIB_INTERP_LINEAR), {2.5f, 0.5f, 0}),
          "Point outside the face projects onto the nearest edge");
    mesh.triangle_count = 0;
    CHECK(near(geomlib_interpolate_vector(&mesh, &vector, {}, GEOMLIB_INTERP_LINEAR), {}),
          "An empty surface returns zero");
    geomlib_scalar_field_free(&scalar);
    geomlib_vector_field_free(&vector);
    geomlib_trimesh_free(&mesh);
    return 0;
}

static int test_tetra_cfd() {
    geomlib_tetramesh mesh;
    geomlib_tetramesh_init(&mesh, 6, 3);
    for (auto p : {geomlib_vec3{0, 0, 0}, geomlib_vec3{2, 0, 0}, geomlib_vec3{0, 3, 0},
                   geomlib_vec3{0, 0, 1}, geomlib_vec3{0, 0, -2}, geomlib_vec3{0, 0, 3}})
        geomlib_tetramesh_add_vertex(&mesh, p);
    geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 2, 3);
    geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 2, 4);
    auto scalar = geomlib_scalar_field_create(mesh.vertex_count);
    for (size_t i = 0; i < scalar.count; i++) {
        auto p = mesh.vertices[i];
        scalar.values[i] = 10 + 2 * p.x - 3 * p.y + 4 * p.z;
    }
    for (auto method : {GEOMLIB_GRADIENT_GREEN_GAUSS, GEOMLIB_GRADIENT_GREEN_GAUSS_NODE,
                       GEOMLIB_GRADIENT_LEAST_SQUARES}) {
        geomlib_vector_field gradients = {};
        CHECK(geomlib_compute_gradient_tetra(&mesh, &scalar, method, &gradients) == GEOMLIB_MESH_SUCCESS &&
              gradients.count == 2, "One gradient should be allocated per cell");
        CHECK(near(gradients.values[0], {2, -3, 4}) && near(gradients.values[1], {2, -3, 4}),
              "An affine gradient is exact for both tetrahedron orientations");
        geomlib_vector_field_free(&gradients);
    }

    geomlib_fvm_mesh fvm = {};
    CHECK(geomlib_fvm_build(&fvm, &mesh) == GEOMLIB_MESH_SUCCESS &&
          fvm.face_count == 7 && fvm.cell_count == 2 && fvm.tetramesh == &mesh,
          "Two tetrahedra should share exactly one face");
    CHECK(near(fvm.cell_volumes[0], 1) && near(fvm.cell_volumes[1], 2) &&
          near(fvm.cell_centers[0], {0.5f, 0.75f, 0.25f}), "Cell centers and volumes");
    unsigned interior = 0;
    for (size_t i = 0; i < fvm.face_count; i++) {
        auto &face = fvm.faces[i];
        auto outward = geomlib_vec3_sub(face.center, fvm.cell_centers[face.cell_left]);
        CHECK(face.vertex_count == 3 && face.area > 0 && near(geomlib_vec3_length(face.normal), 1) &&
              geomlib_vec3_dot(outward, face.normal) > 0, "Face has an outward unit normal and positive area");
        if (face.cell_right != UINT32_MAX) {
            interior++;
            CHECK(face.cell_left == 0 && face.cell_right == 1 && near(face.area, 3) &&
                  fvm.cell_faces[0][3] == i && fvm.cell_faces[1][3] == i,
                  "Shared face and cell-to-face links must agree");
        }
    }
    CHECK(interior == 1 && fvm.cell_face_count[0] == 4 && fvm.cell_face_count[1] == 4,
          "Exactly one interior face and four faces per cell");
    geomlib_fvm_free(&fvm);
    auto thresholds = geomlib_cfd_quality_thresholds_default();
    auto quality = geomlib_compute_cfd_quality(&mesh, 0, &thresholds);
    CHECK(quality.is_valid && near(quality.expansion_ratio, 2) &&
          near(quality.non_orthogonality, 0), "CFD quality measures shared-face geometry");
    thresholds.max_expansion_ratio = 1.5f;
    size_t failures;
    CHECK(!geomlib_check_cfd_mesh_quality(&mesh, &thresholds, &failures) && failures == 2,
          "The configured expansion ratio is enforced");


    geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 2, 5);
    CHECK(geomlib_fvm_build(&fvm, &mesh) == GEOMLIB_MESH_ERROR_INVALID && !fvm.faces,
          "An oversubscribed face must be rejected");
    mesh.tetrahedron_count = 2;
    mesh.tetrahedra[1].vertices[3] = 2;
    geomlib_vector_field gradients = {};
    CHECK(geomlib_compute_gradient_tetra(&mesh, &scalar, GEOMLIB_GRADIENT_GREEN_GAUSS,
          &gradients) == GEOMLIB_MESH_ERROR_DEGENERATE && !gradients.values,
          "Singular cells must leave gradient output unchanged");
    geomlib_scalar_field_free(&scalar);
    geomlib_tetramesh_free(&mesh);
    return 0;
}

static bool read_vtk(const char *path, const geomlib_tetramesh &mesh) {
    std::ifstream in(path);
    std::string token;
    for (int i = 0; i < 4; i++) std::getline(in, token);
    if (token != "DATASET UNSTRUCTURED_GRID") return false;
    size_t count;
    in >> token >> count;
    if (token != "POINTS" || count != mesh.vertex_count) return false;
    in >> token;
    if (token != "float") return false;
    for (size_t i = 0; i < count; i++) {
        geomlib_vec3 p;
        if (!(in >> p.x >> p.y >> p.z) || !near(p, mesh.vertices[i])) return false;
    }
    size_t entries;
    in >> token >> count >> entries;
    if (token != "CELLS" || count != mesh.tetrahedron_count || entries != count * 5) return false;
    for (size_t i = 0; i < count; i++) {
        unsigned arity;
        if (!(in >> arity) || arity != 4) return false;
        for (uint32_t expected : mesh.tetrahedra[i].vertices) {
            uint32_t id;
            if (!(in >> id) || id != expected) return false;
        }
    }
    in >> token >> count;
    if (token != "CELL_TYPES" || count != mesh.tetrahedron_count) return false;
    for (size_t i = 0; i < count; i++) {
        int type;
        if (!(in >> type) || type != 10) return false;
    }
    return true;
}

static bool read_gmsh(const char *path, const geomlib_tetramesh &mesh) {
    std::ifstream in(path);
    std::string token;
    for (int i = 0; i < 4; i++) std::getline(in, token);
    if (token != "$Nodes") return false;
    size_t count;
    if (!(in >> count) || count != mesh.vertex_count) return false;
    for (size_t i = 0; i < count; i++) {
        size_t id;
        geomlib_vec3 p;
        if (!(in >> id >> p.x >> p.y >> p.z) || id != i + 1 || !near(p, mesh.vertices[i])) return false;
    }
    in >> token;
    if (token != "$EndNodes") return false;
    in >> token >> count;
    if (token != "$Elements" || count != mesh.tetrahedron_count) return false;
    for (size_t i = 0; i < count; i++) {
        size_t id;
        int type, tags;
        if (!(in >> id >> type >> tags) || id != i + 1 || type != 4 || tags != 0) return false;
        for (uint32_t expected : mesh.tetrahedra[i].vertices) {
            uint32_t vertex;
            if (!(in >> vertex) || vertex != expected + 1) return false;
        }
    }
    return (in >> token) && token == "$EndElements";
}

static int test_conforming_tetrahedralization() {
    geomlib_trimesh surface;
    geomlib_trimesh_init(&surface, 6, 8);
    for (auto p : {geomlib_vec3{1, 0, 0}, geomlib_vec3{0, 1, 0}, geomlib_vec3{-1, 0, 0},
                   geomlib_vec3{0, -1, 0}, geomlib_vec3{0, 0, 1}, geomlib_vec3{0, 0, -1}})
        geomlib_trimesh_add_vertex(&surface, p, nullptr, nullptr);
    for (uint32_t i = 0; i < 4; i++) {
        geomlib_trimesh_add_triangle(&surface, i, (i + 1) % 4, 4);
        geomlib_trimesh_add_triangle(&surface, (i + 1) % 4, i, 5);
    }

    geomlib_tetramesh mesh = {};
    auto params = geomlib_tetra_params_default(GEOMLIB_TETRA_DELAUNAY);
    CHECK(geomlib_tetrahedralize_surface(&surface, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.tetrahedron_count == 8 && mesh.vertex_count == 7 && mesh.boundary_face_count == 8,
          "An octahedron must produce its own boundary and eight cells");
    CHECK(!mesh.is_boundary_vertex[6], "Fan center is an interior vertex");
    for (size_t i = 0; i < surface.vertex_count; i++) CHECK(mesh.is_boundary_vertex[i], "Mark the boundary rim");

    geomlib_tetramesh_stats stats;
    geomlib_tetramesh_compute_stats(&mesh, &stats);
    CHECK(near(stats.total_volume, 4.0f / 3.0f), "Tetrahedral volume equals the octahedron volume");
    auto quality = geomlib_compute_tetra_quality(&mesh, 0);
    CHECK(near(quality.max_dihedral_angle, 90) && quality.min_dihedral_angle > 54 &&
          quality.min_dihedral_angle < 55, "Report the actual dihedral extrema");
    CHECK(geomlib_reverse_orientation(&surface) == GEOMLIB_MESH_SUCCESS, "Reverse the input shell");
    const geomlib_vec3 points[] = {{0, 0, 0}, {0.2f, 0, 0.3f}};
    CHECK(geomlib_tetrahedralize_constrained(&surface, points, 2, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.vertex_count == 8, "Insert a point on an internal fan face without changing the shell");
    CHECK(geomlib_reverse_orientation(&surface) == GEOMLIB_MESH_SUCCESS, "Restore outward reference winding");
    geomlib_tetramesh_compute_stats(&mesh, &stats);
    CHECK(near(stats.total_volume, 4.0f / 3.0f), "Interior insertion preserves the volume");

    geomlib_trimesh extracted;
    geomlib_trimesh_init(&extracted, 1, 1);
    geomlib_trimesh_add_vertex(&extracted, {-2, 0, 0}, nullptr, nullptr);
    CHECK(geomlib_tetramesh_extract_boundary(&mesh, &extracted) == GEOMLIB_MESH_SUCCESS &&
          extracted.triangle_count == surface.triangle_count, "Extract the conforming shell");
    for (size_t i = 0; i < surface.triangle_count; i++) {
        const uint32_t *expected = &surface.indices[i * 3];
        bool found = false;
        for (size_t k = 0; k < extracted.triangle_count; k++) {
            const uint32_t *actual = &extracted.indices[k * 3];
            for (int rotation = 0; rotation < 3; rotation++) {
                found |= actual[rotation] == expected[0] + 1 &&
                         actual[(rotation + 1) % 3] == expected[1] + 1 &&
                         actual[(rotation + 2) % 3] == expected[2] + 1;
            }
        }
        CHECK(found, "Boundary preserves every original oriented face, including the append offset");
    }
    geomlib_trimesh_free(&extracted);

    const char *vtk = "geomlib_cfd_export_test.vtk";
    const char *msh = "geomlib_cfd_export_test.msh";
    bool round_trip = geomlib_tetramesh_export_vtk(&mesh, vtk) == GEOMLIB_MESH_SUCCESS &&
                      geomlib_tetramesh_export_gmsh(&mesh, msh) == GEOMLIB_MESH_SUCCESS &&
                      read_vtk(vtk, mesh) && read_gmsh(msh, mesh);
    std::remove(vtk);
    std::remove(msh);
    CHECK(round_trip, "Parse emitted VTK and Gmsh coordinates, cell types, and connectivity");

    size_t cells = mesh.tetrahedron_count;
    geomlib_vec3 exterior = {3, 0, 0};
    CHECK(geomlib_tetrahedralize_constrained(&surface, &exterior, 1, &mesh) ==
          GEOMLIB_MESH_ERROR_INVALID && mesh.tetrahedron_count == cells, "Reject outside points atomically");
    surface.triangle_count--;
    CHECK(geomlib_tetrahedralize_surface(&surface, &params, &mesh) == GEOMLIB_MESH_ERROR_INVALID &&
          mesh.tetrahedron_count == cells, "Reject an open surface and preserve the previous mesh");
    geomlib_tetramesh_add_vertex(&mesh, exterior);
    CHECK(!mesh.is_boundary_vertex && !mesh.boundary_faces, "Topology changes invalidate cached tags");
    geomlib_tetramesh_mark_boundary(&mesh);
    CHECK(mesh.boundary_face_count == 8 && !mesh.is_boundary_vertex[mesh.vertex_count - 1],
          "Retag boundary without marking an isolated point");
    geomlib_tetramesh_free(&mesh);
    geomlib_trimesh_free(&surface);
    return 0;
}

static std::vector<std::array<uint32_t, 3>> ordered_faces(const uint32_t *faces, size_t count) {
    std::vector<std::array<uint32_t, 3>> keys(count);
    for (size_t i = 0; i < count; i++) {
        const uint32_t *v = &faces[3 * i];
        size_t first = std::min_element(v, v + 3) - v;
        keys[i] = {v[first], v[(first + 1) % 3], v[(first + 2) % 3]};
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

static double triangle_area(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c) {
    geomlib_vec3 cross = geomlib_vec3_cross(geomlib_vec3_sub(b, a), geomlib_vec3_sub(c, a));
    return 0.5 * geomlib_vec3_length(cross);
}

static double surface_area(const geomlib_trimesh *surface) {
    double area = 0;
    for (size_t i = 0; i < surface->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(surface, i, &a, &b, &c);
        area += triangle_area(a, b, c);
    }
    return area;
}

static double boundary_area(const geomlib_tetramesh *mesh) {
    double area = 0;
    for (size_t i = 0; i < mesh->boundary_face_count; i++) {
        const uint32_t *v = &mesh->boundary_faces[i * 3];
        area += triangle_area(mesh->vertices[v[0]], mesh->vertices[v[1]], mesh->vertices[v[2]]);
    }
    return area;
}

/* Every boundary face must lie inside a single input triangle. */
static bool boundary_subdivides_surface(const geomlib_trimesh *surface,
                                        const geomlib_tetramesh *mesh) {
    const float tolerance = 1e-3f;
    for (size_t f = 0; f < mesh->boundary_face_count; f++) {
        const uint32_t *v = &mesh->boundary_faces[f * 3];
        geomlib_vec3 points[3] = {mesh->vertices[v[0]], mesh->vertices[v[1]], mesh->vertices[v[2]]};
        geomlib_vec3 centroid = geomlib_vec3_mul(
            geomlib_vec3_add(geomlib_vec3_add(points[0], points[1]), points[2]), 1.0f / 3.0f);
        float distance;
        geomlib_find_closest_point_on_mesh(surface, nullptr, centroid, &distance);
        if (distance > tolerance) return false;
        bool contained = false;
        for (size_t t = 0; t < surface->triangle_count && !contained; t++) {
            geomlib_vec3 a, b, c;
            geomlib_trimesh_get_triangle_vertices(surface, t, &a, &b, &c);
            bool inside = true;
            for (int i = 0; i < 3 && inside; i++) {
                geomlib_vec3 nearest = geomlib_closest_point_triangle(points[i], a, b, c);
                inside = geomlib_vec3_length(geomlib_vec3_sub(points[i], nearest)) <= tolerance;
            }
            contained = inside;
        }
        if (!contained) return false;
    }
    return true;
}

static int test_nonstar_shell() {
    geomlib_vec3 lattice[32];
    size_t count = 0;
    for (int z = 0; z < 2; z++)
        for (int y = 0; y < 4; y++)
            for (int x = 0; x < 4; x++)
                lattice[count++] = {float(x), float(y), float(z)};
    geomlib_tetramesh seed = {};
    CHECK(geomlib_tetrahedralize_points(lattice, count, &seed) == GEOMLIB_MESH_SUCCESS,
          "Tetrahedralize a regular lattice");
    size_t out = 0;
    for (size_t i = 0; i < seed.tetrahedron_count; i++) {
        const auto cell = seed.tetrahedra[i];
        const uint32_t *v = cell.vertices;
        auto center = geomlib_tetrahedron_centroid(lattice[v[0]], lattice[v[1]], lattice[v[2]], lattice[v[3]]);
        if (center.y < 1 || center.x < 1 || center.x > 2)
            seed.tetrahedra[out++] = cell;
    }
    seed.tetrahedron_count = out;
    geomlib_trimesh shell = {};
    CHECK(geomlib_tetramesh_extract_boundary(&seed, &shell) == GEOMLIB_MESH_SUCCESS,
          "Extract the U-shaped shell from the lattice");
    geomlib_tetramesh_free(&seed);
    auto params = geomlib_tetra_params_default(GEOMLIB_TETRA_DELAUNAY);
    geomlib_tetramesh filled = {};
    CHECK(geomlib_tetrahedralize_surface(&shell, &params, &filled) == GEOMLIB_MESH_SUCCESS,
          "Fill a connected U-shaped solid whose opposing inner walls rule out a common fan seed");
    CHECK(filled.vertex_count == shell.vertex_count &&
          ordered_faces(filled.boundary_faces, filled.boundary_face_count) ==
          ordered_faces(shell.indices, shell.triangle_count),
          "Fallback must exactly preserve the non-star shell and outward face winding");
    geomlib_tetramesh_stats statistics;
    geomlib_tetramesh_compute_stats(&filled, &statistics);
    CHECK(near(statistics.total_volume, 7), "The fallback should leave the U-shaped notch empty");

    bool changed_diagonal = false;
    for (size_t i = 0; !changed_diagonal && i < shell.triangle_count; i++) {
        uint32_t *first = &shell.indices[i * 3];
        if (*std::max_element(first, first + 3) >= 16) continue;
        for (size_t j = i + 1; !changed_diagonal && j < shell.triangle_count; j++) {
            uint32_t *second = &shell.indices[j * 3];
            if (*std::max_element(second, second + 3) >= 16) continue;
            for (int k = 0; !changed_diagonal && k < 3; k++) {
                uint32_t a = first[k], b = first[(k + 1) % 3];
                if (std::abs(int(a % 4) - int(b % 4)) != 1 ||
                    std::abs(int(a / 4) - int(b / 4)) != 1) continue;
                for (int l = 0; l < 3; l++) {
                    if (second[l] != b || second[(l + 1) % 3] != a) continue;
                    uint32_t c = first[(k + 2) % 3], d = second[(l + 2) % 3];
                    first[0] = c; first[1] = a; first[2] = d;
                    second[0] = c; second[1] = d; second[2] = b;
                    changed_diagonal = true;
                    break;
                }
            }
        }
    }
    CHECK(changed_diagonal, "The fixture has a square boundary patch to retriangulate");
    size_t old_count = filled.tetrahedron_count;
    CHECK(geomlib_tetrahedralize_surface(&shell, &params, &filled) == GEOMLIB_MESH_SUCCESS,
          "Facet recovery should force the constrained boundary diagonal");
    CHECK(filled.tetrahedron_count != old_count || filled.boundary_face_count > shell.triangle_count,
          "Recovery should insert Steiner points rather than keep the wrong diagonal");
    geomlib_tetramesh_compute_stats(&filled, &statistics);
    CHECK(near(statistics.total_volume, 7), "Recovered shell keeps the U-shaped volume");
    CHECK(filled.boundary_face_count >= shell.triangle_count &&
          std::abs(boundary_area(&filled) - surface_area(&shell)) < 1e-3 * surface_area(&shell),
          "Recovered boundary area matches the input surface");
    CHECK(boundary_subdivides_surface(&shell, &filled),
          "Every recovered boundary face lies inside an input triangle");

    const geomlib_vec3 flat[] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}};
    size_t cells = filled.tetrahedron_count;
    CHECK(geomlib_tetrahedralize_points(flat, 4, &filled) == GEOMLIB_MESH_ERROR_DEGENERATE &&
          filled.tetrahedron_count == cells, "Reject coplanar points without replacing the old result");
    CHECK(geomlib_tetrahedralize_box({0,0,0}, {1,1,1}, &params, &filled) == GEOMLIB_MESH_SUCCESS,
          "Box meshing uses the validated Delaunay result");
    geomlib_tetramesh_compute_stats(&filled, &statistics);
    CHECK(filled.boundary_face_count == 12 && near(statistics.total_volume, 1),
          "A unit box has twelve boundary triangles and unit volume");
    geomlib_tetramesh_free(&filled);
    geomlib_trimesh_free(&shell);
    return 0;
}

static int test_tetrahedral_edits() {
    geomlib_tetramesh mesh = {};
    for (auto p : {geomlib_vec3{1, 0, 0}, geomlib_vec3{-0.5f, 0.8660254f, 0},
                   geomlib_vec3{-0.5f, -0.8660254f, 0}, geomlib_vec3{0, 0, -0.001f},
                   geomlib_vec3{0, 0, 1}})
        CHECK(geomlib_tetramesh_add_vertex(&mesh, p) != UINT32_MAX, "Add fixture vertices");
    geomlib_tetramesh_add_tetrahedron(&mesh, 0, 2, 1, 3);
    geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 2, 4);
    geomlib_tetramesh_mark_boundary(&mesh);
    auto original_faces = ordered_faces(mesh.boundary_faces, mesh.boundary_face_count);
    geomlib_tetramesh_stats before, after;
    geomlib_tetramesh_compute_stats(&mesh, &before);
    CHECK(before.sliver_count == 1, "The bipyramid contains a sliver");
    geomlib_tetramesh optimized = {};
    for (size_t i = 0; i < mesh.vertex_count; i++)
        geomlib_tetramesh_add_vertex(&optimized, mesh.vertices[i]);
    geomlib_tetramesh_add_tetrahedron(&optimized, 0, 2, 1, 3);
    geomlib_tetramesh_add_tetrahedron(&optimized, 0, 1, 2, 4);
    CHECK(geomlib_tetramesh_optimize(&optimized, 1) == GEOMLIB_MESH_SUCCESS &&
          optimized.tetrahedron_count == 3 &&
          original_faces == ordered_faces(optimized.boundary_faces, optimized.boundary_face_count),
          "Optimization uses its flip budget without changing the shell");
    geomlib_tetramesh_free(&optimized);
    CHECK(geomlib_tetramesh_remove_slivers(&mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.tetrahedron_count == 3 && mesh.vertex_count == 5, "A 2-to-3 flip removes the sliver");
    geomlib_tetramesh_compute_stats(&mesh, &after);
    CHECK(after.sliver_count == 0 && near(after.total_volume, before.total_volume) &&
          original_faces == ordered_faces(mesh.boundary_faces, mesh.boundary_face_count),
          "Removing the sliver preserves volume and every boundary face");
    CHECK(geomlib_tetramesh_remove_slivers(&mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.tetrahedron_count == 3, "A second sliver pass has no changes");

    auto params = geomlib_tetra_params_default(GEOMLIB_TETRA_DELAUNAY);
    params.max_edge_length = 0.9f;
    params.max_insertions = 1;
    CHECK(geomlib_tetramesh_refine(&mesh, &params) == GEOMLIB_MESH_SUCCESS &&
          mesh.vertex_count == 6 && mesh.tetrahedron_count == 6 && !mesh.is_boundary_vertex[5],
          "Refine all three cells sharing the long interior edge");
    geomlib_tetramesh_compute_stats(&mesh, &after);
    CHECK(near(after.total_volume, before.total_volume) &&
          original_faces == ordered_faces(mesh.boundary_faces, mesh.boundary_face_count),
          "Interior edge bisection preserves the shell and volume");

    params.preserve_boundary = false;
    params.refine_boundary = true;
    CHECK(geomlib_tetramesh_refine(&mesh, &params) == GEOMLIB_MESH_SUCCESS &&
          mesh.vertex_count == 7 && mesh.is_boundary_vertex[6] &&
          mesh.boundary_face_count == original_faces.size() + 2,
          "Explicit boundary refinement bisects the boundary faces incident to the selected edge");
    geomlib_tetramesh_compute_stats(&mesh, &after);
    CHECK(near(after.total_volume, before.total_volume), "Boundary subdivision preserves volume");

    size_t cells = mesh.tetrahedron_count;
    params.max_edge_length = NAN;
    CHECK(geomlib_tetramesh_refine(&mesh, &params) == GEOMLIB_MESH_ERROR_INVALID &&
          mesh.tetrahedron_count == cells, "Invalid size leaves the input intact");
    CHECK(geomlib_tetramesh_optimize(&mesh, -1) == GEOMLIB_MESH_ERROR_INVALID &&
          geomlib_tetramesh_optimize(&mesh, 0) == GEOMLIB_MESH_SUCCESS,
          "Optimization validates its budget");
    geomlib_tetramesh_free(&mesh);
    return 0;
}

static int test_tetrahedral_quality() {
    geomlib_tetramesh mesh = {};
    const geomlib_vec3 points[] = {{0, 0, 0}, {0, 0, 2}, {1, 0, 1}, {-1, 0.001f, 1}, {0, -1, 1}};
    for (geomlib_vec3 point : points)
        CHECK(geomlib_tetramesh_add_vertex(&mesh, point) != UINT32_MAX, "Add 3-to-2 fixture vertices");
    CHECK(geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 2, 3) != UINT32_MAX &&
          geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 3, 4) != UINT32_MAX &&
          geomlib_tetramesh_add_tetrahedron(&mesh, 0, 1, 4, 2) != UINT32_MAX,
          "Three cells share the 0-1 edge");
    geomlib_tetramesh_mark_boundary(&mesh);
    auto original_faces = ordered_faces(mesh.boundary_faces, mesh.boundary_face_count);
    geomlib_tetramesh_stats before, after;
    geomlib_tetramesh_compute_stats(&mesh, &before);
    CHECK(before.sliver_count > 0, "The fixture contains slivers");

    CHECK(geomlib_tetramesh_remove_slivers(&mesh) == GEOMLIB_MESH_SUCCESS,
          "Sliver removal should succeed");
    geomlib_tetramesh_compute_stats(&mesh, &after);
    CHECK(mesh.tetrahedron_count == 2, "A 3-to-2 flip replaces three cells with two");
    CHECK(near(after.total_volume, before.total_volume) &&
          original_faces == ordered_faces(mesh.boundary_faces, mesh.boundary_face_count),
          "The 3-to-2 flip preserves volume and the boundary face set");
    CHECK(after.sliver_count < before.sliver_count, "The flip removes slivers");
    for (size_t i = 0; i < mesh.tetrahedron_count; i++)
        CHECK(geomlib_compute_tetra_quality(&mesh, (uint32_t)i).volume > 0, "Cells stay positive");
    geomlib_tetramesh_free(&mesh);
    return 0;
}

static int test_tetrahedral_smoothing() {
    geomlib_trimesh surface;
    geomlib_trimesh_init(&surface, 6, 8);
    for (geomlib_vec3 point : {geomlib_vec3{1, 0, 0}, geomlib_vec3{0, 1, 0},
                               geomlib_vec3{-1, 0, 0}, geomlib_vec3{0, -1, 0},
                               geomlib_vec3{0, 0, 1}, geomlib_vec3{0, 0, -1}})
        geomlib_trimesh_add_vertex(&surface, point, nullptr, nullptr);
    for (uint32_t i = 0; i < 4; i++) {
        geomlib_trimesh_add_triangle(&surface, i, (i + 1) % 4, 4);
        geomlib_trimesh_add_triangle(&surface, (i + 1) % 4, i, 5);
    }
    auto params = geomlib_tetra_params_default(GEOMLIB_TETRA_DELAUNAY);
    geomlib_tetramesh mesh = {};
    CHECK(geomlib_tetrahedralize_surface(&surface, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.tetrahedron_count == 8, "Octahedron fixture builds eight cells");
    CHECK(mesh.vertex_count == 7 && !mesh.is_boundary_vertex[6], "The fan center is interior");
    mesh.vertices[6] = {0.5f, 0.2f, 0.1f};
    geomlib_tetramesh_stats before, after;
    geomlib_tetramesh_compute_stats(&mesh, &before);
    auto original_faces = ordered_faces(mesh.boundary_faces, mesh.boundary_face_count);

    CHECK(geomlib_tetramesh_smooth(&mesh, 20) == GEOMLIB_MESH_SUCCESS,
          "Interior smoothing should succeed");
    geomlib_tetramesh_compute_stats(&mesh, &after);
    CHECK(near(after.total_volume, before.total_volume),
          "Smoothing preserves the total volume");
    CHECK(original_faces == ordered_faces(mesh.boundary_faces, mesh.boundary_face_count),
          "Smoothing leaves the boundary faces fixed");
    CHECK(!mesh.is_boundary_vertex[6], "The smoothed vertex stays interior");
    CHECK(geomlib_vec3_length(mesh.vertices[6]) < 0.2f,
          "The center relaxes toward the octahedron centroid");
    CHECK(after.min_quality > before.min_quality, "Smoothing improves the worst cell quality");
    CHECK(geomlib_tetramesh_smooth(&mesh, 0) == GEOMLIB_MESH_SUCCESS &&
          geomlib_tetramesh_smooth(&mesh, -1) == GEOMLIB_MESH_ERROR_INVALID,
          "Smoothing validates its iteration budget");
    geomlib_tetramesh_free(&mesh);
    geomlib_trimesh_free(&surface);
    return 0;
}

static int test_streamline_interpolation() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 2);
    geomlib_trimesh_add_vertex(&mesh, {0, 0, 0}, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, {1, 0, 0}, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, {1, 1, 0}, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, {0, 1, 0}, nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 2, 3);

    auto velocity = geomlib_vector_field_create(mesh.vertex_count);
    for (size_t i = 0; i < mesh.vertex_count; i++)
        velocity.values[i] = {mesh.vertices[i].x, mesh.vertices[i].y, 0};

    geomlib_streamline_params params;
    params.method = GEOMLIB_INTEGRATE_EULER;
    params.step_size = 0.1f;
    params.max_steps = 2;
    params.max_length = 10.0f;
    params.bounds = {{-10, -10, -10}, {10, 10, 10}};

    geomlib_streamline line = geomlib_streamline_create(8);
    CHECK(geomlib_trace_streamline(&mesh, &velocity, {0.4f, 0.2f, 0}, &params, &line) ==
          GEOMLIB_MESH_SUCCESS && line.point_count >= 2,
          "Streamline tracing should produce points");
    /* The velocity field is v = (x, y); the first Euler step must use the
     * barycentric value (0.4, 0.2), not a vertex value. */
    CHECK(near(line.points[1], {0.44f, 0.22f, 0}),
          "Streamline velocity is interpolated barycentrically");

    geomlib_streamline_free(&line);
    free(velocity.values);
    geomlib_trimesh_free(&mesh);
    return 0;
}

int main() {
    const int failures = test_surface_operators() + test_interpolation() + test_tetra_cfd() +
                         test_conforming_tetrahedralization() + test_nonstar_shell() +
                         test_tetrahedral_edits() + test_tetrahedral_quality() +
                         test_tetrahedral_smoothing() + test_streamline_interpolation();
    if (!failures) std::puts("All CFD tests passed.");
    return failures;
}
