/*
 * geomlib - Tetrahedral Mesh Implementation
 * Delaunay tetrahedralization for CFD simulations
 */

#include "geomlib/geomlib_tetramesh.h"
#include "geomlib/geomlib_cfd.h"
#include "geomlib/geomlib_sdf.h"
#include "geomlib/geomlib_montecarlo.h"
#include "geomlib/geomlib_statistics.h"
#include "geomlib/geomlib_topology.h"
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Basic Tetrahedral Mesh Operations
 * ============================================================================ */

void geomlib_tetramesh_init(
    geomlib_tetramesh *mesh,
    size_t vertex_capacity,
    size_t tetrahedron_capacity
) {
    if (!mesh) return;
    *mesh = (geomlib_tetramesh){0};
    if (vertex_capacity <= SIZE_MAX / sizeof(*mesh->vertices) && vertex_capacity <= UINT32_MAX)
        mesh->vertices = malloc(vertex_capacity * sizeof(*mesh->vertices));
    if (tetrahedron_capacity <= SIZE_MAX / sizeof(*mesh->tetrahedra) && tetrahedron_capacity <= UINT32_MAX)
        mesh->tetrahedra = malloc(tetrahedron_capacity * sizeof(*mesh->tetrahedra));
    if (mesh->vertices) mesh->vertex_capacity = vertex_capacity;
    if (mesh->tetrahedra) mesh->tetrahedron_capacity = tetrahedron_capacity;
}

void geomlib_tetramesh_free(geomlib_tetramesh *mesh) {
    if (!mesh) return;

    free(mesh->vertices);
    free(mesh->tetrahedra);
    free(mesh->boundary_faces);
    free(mesh->is_boundary_vertex);

    memset(mesh, 0, sizeof(geomlib_tetramesh));
}

static void invalidate_boundary(geomlib_tetramesh *mesh) {
    free(mesh->boundary_faces);
    free(mesh->is_boundary_vertex);
    mesh->boundary_faces = NULL;
    mesh->is_boundary_vertex = NULL;
    mesh->boundary_face_count = 0;
}

uint32_t geomlib_tetramesh_add_vertex(
    geomlib_tetramesh *mesh,
    geomlib_vec3 position
) {
    if (!mesh || mesh->vertex_count >= UINT32_MAX) return UINT32_MAX;
    if (mesh->vertex_count == mesh->vertex_capacity) {
        size_t limit = SIZE_MAX / sizeof(*mesh->vertices);
        if (limit > UINT32_MAX) limit = UINT32_MAX;
        if (mesh->vertex_capacity >= limit) return UINT32_MAX;
        size_t next = mesh->vertex_capacity > limit / 2 ? limit : mesh->vertex_capacity * 2;
        if (!next) next = 8;
        geomlib_vec3 *vertices = realloc(mesh->vertices, next * sizeof(*vertices));
        if (!vertices) return UINT32_MAX;
        mesh->vertices = vertices;
        mesh->vertex_capacity = next;
    }

    invalidate_boundary(mesh);
    mesh->vertices[mesh->vertex_count] = position;
    return mesh->vertex_count++;
}

uint32_t geomlib_tetramesh_add_tetrahedron(
    geomlib_tetramesh *mesh,
    uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3
) {
    if (!mesh || mesh->tetrahedron_count >= UINT32_MAX) return UINT32_MAX;
    if (mesh->tetrahedron_count == mesh->tetrahedron_capacity) {
        size_t limit = SIZE_MAX / sizeof(*mesh->tetrahedra);
        if (limit > UINT32_MAX) limit = UINT32_MAX;
        if (mesh->tetrahedron_capacity >= limit) return UINT32_MAX;
        size_t next = mesh->tetrahedron_capacity > limit / 2 ? limit : mesh->tetrahedron_capacity * 2;
        if (!next) next = 8;
        geomlib_tetrahedron *cells = realloc(mesh->tetrahedra, next * sizeof(*cells));
        if (!cells) return UINT32_MAX;
        mesh->tetrahedra = cells;
        mesh->tetrahedron_capacity = next;
    }

    geomlib_tetrahedron tet;
    tet.vertices[0] = v0;
    tet.vertices[1] = v1;
    tet.vertices[2] = v2;
    tet.vertices[3] = v3;

    invalidate_boundary(mesh);
    mesh->tetrahedra[mesh->tetrahedron_count] = tet;
    return mesh->tetrahedron_count++;
}

/* Default parameters */
geomlib_tetra_params geomlib_tetra_params_default(geomlib_tetra_method method) {
    geomlib_tetra_params params;
    params.method = method;
    params.max_edge_length = 0.0f;        /* Unlimited */
    params.min_edge_length = 0.001f;
    params.element_size = 1.0f;

    params.min_dihedral_angle = 10.0f;    /* Degrees */
    params.max_aspect_ratio = 10.0f;
    params.optimize_quality = true;

    params.preserve_boundary = true;
    params.refine_boundary = false;

    params.max_iterations = 100;
    params.max_insertions = 10000;

    return params;
}

/* ============================================================================
 * Geometric Predicates and Utilities
 * ============================================================================ */

/* Compute signed volume of tetrahedron (1/6 of scalar triple product) */
float geomlib_tetrahedron_volume(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
) {
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 e3 = geomlib_vec3_sub(v3, v0);

    /* Volume = (1/6) * |e1 · (e2 × e3)| */
    geomlib_vec3 cross = geomlib_vec3_cross(e2, e3);
    float dot = geomlib_vec3_dot(e1, cross);
    return fabsf(dot) / 6.0f;
}

/* Compute circumsphere of tetrahedron */
void geomlib_tetrahedron_circumsphere(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3,
    geomlib_vec3 *center,
    float *radius
) {
    /* Using barycentric approach */
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 e3 = geomlib_vec3_sub(v3, v0);

    float l1 = geomlib_vec3_dot(e1, e1);
    float l2 = geomlib_vec3_dot(e2, e2);
    float l3 = geomlib_vec3_dot(e3, e3);

    geomlib_vec3 c1 = geomlib_vec3_cross(e2, e3);
    geomlib_vec3 c2 = geomlib_vec3_cross(e3, e1);
    geomlib_vec3 c3 = geomlib_vec3_cross(e1, e2);

    float det = 2.0f * geomlib_vec3_dot(e1, c1);

    if (fabsf(det) < 1e-10f) {
        /* Degenerate tetrahedron */
        *center = v0;
        *radius = FLT_MAX;
        return;
    }

    geomlib_vec3 offset = geomlib_vec3_create(
        (l1 * c1.x + l2 * c2.x + l3 * c3.x) / det,
        (l1 * c1.y + l2 * c2.y + l3 * c3.y) / det,
        (l1 * c1.z + l2 * c2.z + l3 * c3.z) / det
    );

    *center = geomlib_vec3_add(v0, offset);
    *radius = geomlib_vec3_length(offset);
}

static geomlib_vec3d tetra_point(geomlib_vec3 point) {
    return (geomlib_vec3d){point.x, point.y, point.z};
}

static geomlib_vec3 tetra_point_float(geomlib_vec3d point) {
    return (geomlib_vec3){(float)point.x, (float)point.y, (float)point.z};
}

static bool finite_point(geomlib_vec3 point) {
    return isfinite(point.x) && isfinite(point.y) && isfinite(point.z);
}

static double oriented_volume6(geomlib_vec3d a, geomlib_vec3d b, geomlib_vec3d c, geomlib_vec3d d) {
    return geomlib_vec3d_dot(geomlib_vec3d_sub(b, a),
        geomlib_vec3d_cross(geomlib_vec3d_sub(c, a), geomlib_vec3d_sub(d, a)));
}

/* Check orientation (sign of volume) */
static float orientation(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
) {
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 e3 = geomlib_vec3_sub(v3, v0);

    geomlib_vec3 cross = geomlib_vec3_cross(e2, e3);
    return geomlib_vec3_dot(e1, cross);
}

/* Check if point is inside tetrahedron using barycentric coordinates */
bool geomlib_point_in_tetrahedron(
    geomlib_vec3 point,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
) {
    /* Check if point is on the same side as each vertex */
    float o0 = orientation(v1, v2, v3, point);
    float o1 = orientation(v0, v3, v2, point);
    float o2 = orientation(v0, v1, v3, point);
    float o3 = orientation(v0, v2, v1, point);

    /* Check if all have same sign (including zero for boundary) */
    bool all_pos = (o0 >= 0 && o1 >= 0 && o2 >= 0 && o3 >= 0);
    bool all_neg = (o0 <= 0 && o1 <= 0 && o2 <= 0 && o3 <= 0);

    return all_pos || all_neg;
}

/* ============================================================================
 * Bowyer-Watson Delaunay Tetrahedralization
 * ============================================================================ */

static int tag_boundary(geomlib_tetramesh *mesh);
static double cell_volume6(const geomlib_tetramesh *mesh, const geomlib_tetrahedron *cell);

typedef struct tetra_face_key {
    uint32_t v[3];
} tetra_face_key;

static int compare_face_keys(const void *left, const void *right) {
    const tetra_face_key *a = left, *b = right;
    for (int i = 0; i < 3; i++) {
        if (a->v[i] != b->v[i]) return a->v[i] < b->v[i] ? -1 : 1;
    }
    return 0;
}

static tetra_face_key unordered_face(uint32_t a, uint32_t b, uint32_t c) {
    if (a > b) { uint32_t temp = a; a = b; b = temp; }
    if (b > c) { uint32_t temp = b; b = c; c = temp; }
    if (a > b) { uint32_t temp = a; a = b; b = temp; }
    return (tetra_face_key){{a, b, c}};
}

static bool sphere_contains(geomlib_vec3d point, const geomlib_vec3d *positions,
                            const geomlib_tetrahedron *tetra) {
    geomlib_vec3d origin = positions[tetra->vertices[0]];
    geomlib_vec3d edges[3];
    for (int i = 0; i < 3; i++) edges[i] = geomlib_vec3d_sub(positions[tetra->vertices[i + 1]], origin);
    geomlib_vec3d cofactor[3] = {
        geomlib_vec3d_cross(edges[1], edges[2]), geomlib_vec3d_cross(edges[2], edges[0]),
        geomlib_vec3d_cross(edges[0], edges[1])
    };
    double determinant = geomlib_vec3d_dot(edges[0], cofactor[0]);
    if (!determinant) return false;
    geomlib_vec3d center = {0};
    for (int i = 0; i < 3; i++) {
        center = geomlib_vec3d_add(center,
            geomlib_vec3d_mul(cofactor[i], geomlib_vec3d_dot(edges[i], edges[i])));
    }
    center = geomlib_vec3d_mul(center, 0.5 / determinant);
    geomlib_vec3d delta = geomlib_vec3d_sub(geomlib_vec3d_sub(point, origin), center);
    double distance2 = geomlib_vec3d_dot(delta, delta), radius2 = geomlib_vec3d_dot(center, center);
    return distance2 <= radius2 + 32 * DBL_EPSILON * fmax(distance2, radius2);
}

static int insert_delaunay_point(geomlib_tetramesh *mesh, const geomlib_vec3d *positions,
                                 uint32_t index) {
    size_t cells = mesh->tetrahedron_count;
    if (cells > SIZE_MAX / (4 * sizeof(tetra_face_key))) return GEOMLIB_MESH_ERROR_MEMORY;
    tetra_face_key *faces = malloc(cells * 4 * sizeof(*faces));
    if (!faces) return GEOMLIB_MESH_ERROR_MEMORY;

    size_t kept = 0, count = 0;
    for (size_t i = 0; i < cells; i++) {
        geomlib_tetrahedron tetra = mesh->tetrahedra[i];
        if (!sphere_contains(positions[index], positions, &tetra)) {
            mesh->tetrahedra[kept++] = tetra;
            continue;
        }
        const uint32_t *v = tetra.vertices;
        for (int j = 0; j < 4; j++)
            faces[count++] = unordered_face(v[(j + 1) % 4], v[(j + 2) % 4], v[(j + 3) % 4]);
    }
    mesh->tetrahedron_count = kept;
    qsort(faces, count, sizeof(*faces), compare_face_keys);
    int result = count ? GEOMLIB_MESH_SUCCESS : GEOMLIB_MESH_ERROR_DEGENERATE;
    for (size_t start = 0; start < count;) {
        size_t end = start + 1;
        while (end < count && !compare_face_keys(&faces[start], &faces[end])) end++;
        if (end - start > 2) {
            result = GEOMLIB_MESH_ERROR_DEGENERATE;
            break;
        }
        if (end - start == 1) {
            uint32_t a = faces[start].v[0], b = faces[start].v[1], c = faces[start].v[2];
            double volume = oriented_volume6(positions[a], positions[b], positions[c], positions[index]);
            if (volume && geomlib_tetramesh_add_tetrahedron(mesh, a, volume > 0 ? b : c,
                    volume > 0 ? c : b, index) == UINT32_MAX) {
                result = GEOMLIB_MESH_ERROR_MEMORY;
                break;
            }
        }
        start = end;
    }
    free(faces);
    return result;
}

int geomlib_tetrahedralize_points(const geomlib_vec3 *points, size_t count,
                                  geomlib_tetramesh *tetramesh) {
    if (!points || !tetramesh || count < 4 || count > UINT32_MAX - 4 ||
        count > SIZE_MAX / sizeof(geomlib_vec3d) - 4) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_vec3 min = points[0], max = points[0];
    for (size_t i = 0; i < count; i++) {
        geomlib_vec3 p = points[i];
        if (!finite_point(p)) return GEOMLIB_MESH_ERROR_INVALID;
        min = (geomlib_vec3){fminf(min.x, p.x), fminf(min.y, p.y), fminf(min.z, p.z)};
        max = (geomlib_vec3){fmaxf(max.x, p.x), fmaxf(max.y, p.y), fmaxf(max.z, p.z)};
    }
    geomlib_vec3d extent = geomlib_vec3d_sub(tetra_point(max), tetra_point(min));
    double scale = fmax(extent.x, fmax(extent.y, extent.z));
    if (!scale) return GEOMLIB_MESH_ERROR_DEGENERATE;

    geomlib_tetramesh mesh;
    geomlib_tetramesh_init(&mesh, count, count);
    geomlib_vec3d *positions = malloc((count + 4) * sizeof(*positions));
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (!mesh.vertices || !mesh.tetrahedra || !positions) goto finish;

    memcpy(mesh.vertices, points, count * sizeof(*points));
    mesh.vertex_count = count;
    for (size_t i = 0; i < count; i++)
        positions[i] = geomlib_vec3d_mul(geomlib_vec3d_sub(tetra_point(points[i]), tetra_point(min)), 1 / scale);
    uint32_t first = (uint32_t)count;
    positions[first] = (geomlib_vec3d){-16, -16, -16};
    positions[first + 1] = (geomlib_vec3d){64, -16, -16};
    positions[first + 2] = (geomlib_vec3d){-16, 64, -16};
    positions[first + 3] = (geomlib_vec3d){-16, -16, 64};
    if (geomlib_tetramesh_add_tetrahedron(&mesh, first, first + 1, first + 2, first + 3) == UINT32_MAX)
        goto finish;

    for (size_t i = 0; i < count; i++) {
        bool duplicate = false;
        for (size_t j = 0; j < i && !duplicate; j++)
            duplicate = points[i].x == points[j].x && points[i].y == points[j].y && points[i].z == points[j].z;
        if (duplicate) continue;
        result = insert_delaunay_point(&mesh, positions, (uint32_t)i);
        if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    }

    size_t kept = 0;
    for (size_t i = 0; i < mesh.tetrahedron_count; i++) {
        geomlib_tetrahedron cell = mesh.tetrahedra[i];
        if (cell.vertices[0] < first && cell.vertices[1] < first &&
            cell.vertices[2] < first && cell.vertices[3] < first)
            mesh.tetrahedra[kept++] = cell;
    }
    mesh.tetrahedron_count = kept;
    result = kept ? tag_boundary(&mesh) : GEOMLIB_MESH_ERROR_DEGENERATE;

finish:
    free(positions);
    if (result == GEOMLIB_MESH_SUCCESS) {
        geomlib_tetramesh_free(tetramesh);
        *tetramesh = mesh;
    } else {
        geomlib_tetramesh_free(&mesh);
    }
    return result;
}

/* Generate tetrahedral mesh from bounding box */
int geomlib_tetrahedralize_box(
    geomlib_vec3 min_corner,
    geomlib_vec3 max_corner,
    const geomlib_tetra_params *params,
    geomlib_tetramesh *tetramesh
) {
    if (!params || !tetramesh || !finite_point(min_corner) || !finite_point(max_corner) ||
        !(params->element_size > 0) || !isfinite(params->element_size)) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_vec3d extent = geomlib_vec3d_sub(tetra_point(max_corner), tetra_point(min_corner));
    const double lengths[] = {extent.x, extent.y, extent.z};
    size_t count = 1, resolution[3];
    for (int i = 0; i < 3; i++) {
        if (!(lengths[i] > 0)) return GEOMLIB_MESH_ERROR_DEGENERATE;
        double n = fmax(2, floor(lengths[i] / params->element_size) + 1);
        if (n > UINT32_MAX - 4) return GEOMLIB_MESH_ERROR_MEMORY;
        resolution[i] = (size_t)n;
        if (count > SIZE_MAX / sizeof(geomlib_vec3) / resolution[i] ||
            count > (UINT32_MAX - 4) / resolution[i]) return GEOMLIB_MESH_ERROR_MEMORY;
        count *= resolution[i];
    }

    geomlib_vec3 *points = malloc(count * sizeof(*points));
    if (!points) return GEOMLIB_MESH_ERROR_MEMORY;
    size_t index = 0;
    for (size_t z = 0; z < resolution[2]; z++) {
        for (size_t y = 0; y < resolution[1]; y++) {
            for (size_t x = 0; x < resolution[0]; x++) {
                points[index++] = (geomlib_vec3){
                    (float)(min_corner.x + extent.x * (x / (double)(resolution[0] - 1))),
                    (float)(min_corner.y + extent.y * (y / (double)(resolution[1] - 1))),
                    (float)(min_corner.z + extent.z * (z / (double)(resolution[2] - 1)))
                };
            }
        }
    }

    /* Tetrahedralize */
    int result = geomlib_tetrahedralize_points(points, count, tetramesh);

    free(points);
    return result;
}

static bool supports_surface(const geomlib_trimesh *surface, double winding, geomlib_vec3 seed,
                             double tolerance) {
    if (!finite_point(seed)) return false;
    geomlib_vec3d point = tetra_point(seed);
    for (size_t i = 0; i < surface->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(surface, i, &a, &b, &c);
        geomlib_vec3d origin = tetra_point(a);
        geomlib_vec3d normal = geomlib_vec3d_cross(
            geomlib_vec3d_sub(tetra_point(b), origin), geomlib_vec3d_sub(tetra_point(c), origin));
        double length = geomlib_vec3d_length(normal);
        double height = winding * geomlib_vec3d_dot(normal, geomlib_vec3d_sub(point, origin));
        if (!(height < -tolerance * length)) return false;
    }
    return true;
}

static bool cell_coordinates(const geomlib_tetramesh *mesh, size_t cell, geomlib_vec3 points[4]) {
    if (!mesh || !mesh->vertices || !mesh->tetrahedra || cell >= mesh->tetrahedron_count)
        return false;
    const uint32_t *indices = mesh->tetrahedra[cell].vertices;
    for (int i = 0; i < 4; i++) {
        if (indices[i] >= mesh->vertex_count || !finite_point(mesh->vertices[indices[i]])) return false;
        points[i] = mesh->vertices[indices[i]];
    }
    return true;
}

/* Split every containing cell, including both sides of a shared facet. */
static int insert_interior_vertex(geomlib_tetramesh *mesh, geomlib_vec3 point) {
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        geomlib_vec3 p = mesh->vertices[i];
        if (point.x == p.x && point.y == p.y && point.z == p.z) return GEOMLIB_MESH_ERROR_INVALID;
    }
    uint32_t vertex = geomlib_tetramesh_add_vertex(mesh, point);
    if (vertex == UINT32_MAX) return GEOMLIB_MESH_ERROR_MEMORY;
    geomlib_vec3d p = tetra_point(point);
    bool contained = false;
    size_t cells = mesh->tetrahedron_count;
    for (size_t cell = 0; cell < cells; cell++) {
        geomlib_vec3 points[4];
        if (!cell_coordinates(mesh, cell, points)) return GEOMLIB_MESH_ERROR_INVALID;
        geomlib_vec3d corners[4];
        for (int i = 0; i < 4; i++) corners[i] = tetra_point(points[i]);
        double volume = oriented_volume6(corners[0], corners[1], corners[2], corners[3]);
        double weights[4] = {
            oriented_volume6(p, corners[1], corners[2], corners[3]) / volume,
            oriented_volume6(corners[0], p, corners[2], corners[3]) / volume,
            oriented_volume6(corners[0], corners[1], p, corners[3]) / volume,
            oriented_volume6(corners[0], corners[1], corners[2], p) / volume
        };
        const double tolerance = 64 * DBL_EPSILON;
        if (weights[0] < -tolerance || weights[1] < -tolerance ||
            weights[2] < -tolerance || weights[3] < -tolerance) continue;
        bool replaced = false;
        geomlib_tetrahedron parent = mesh->tetrahedra[cell];
        for (int i = 0; i < 4; i++) {
            if (weights[i] <= tolerance) continue;
            geomlib_tetrahedron child = parent;
            child.vertices[i] = vertex;
            if (!replaced) mesh->tetrahedra[cell] = child;
            else if (geomlib_tetramesh_add_tetrahedron(mesh, child.vertices[0], child.vertices[1],
                child.vertices[2], child.vertices[3]) == UINT32_MAX) return GEOMLIB_MESH_ERROR_MEMORY;
            replaced = true;
        }
        contained |= replaced;
    }
    return contained ? GEOMLIB_MESH_SUCCESS : GEOMLIB_MESH_ERROR_INVALID;
}

static tetra_face_key oriented_face(const uint32_t face[3], bool reverse) {
    int first = face[1] < face[0] ? 1 : 0;
    if (face[2] < face[first]) first = 2;
    return (tetra_face_key){{face[first], face[(first + (reverse ? 2 : 1)) % 3],
                             face[(first + (reverse ? 1 : 2)) % 3]}};
}

static int verify_surface_boundary(const geomlib_trimesh *surface, bool reverse,
                                   const geomlib_tetramesh *mesh) {
    size_t count = surface->triangle_count;
    if (mesh->boundary_face_count != count) return GEOMLIB_MESH_ERROR_UNSUPPORTED;
    if (count > SIZE_MAX / sizeof(tetra_face_key)) return GEOMLIB_MESH_ERROR_MEMORY;
    tetra_face_key *expected = malloc(count * sizeof(*expected));
    tetra_face_key *actual = malloc(count * sizeof(*actual));
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (expected && actual) {
        for (size_t i = 0; i < count; i++) {
            expected[i] = oriented_face(&surface->indices[i * 3], reverse);
            actual[i] = oriented_face(&mesh->boundary_faces[i * 3], false);
        }
        qsort(expected, count, sizeof(*expected), compare_face_keys);
        qsort(actual, count, sizeof(*actual), compare_face_keys);
        result = GEOMLIB_MESH_SUCCESS;
        for (size_t i = 0; i < count; i++) {
            if (!compare_face_keys(&expected[i], &actual[i])) continue;
            result = GEOMLIB_MESH_ERROR_UNSUPPORTED;
            break;
        }
    }
    free(expected);
    free(actual);
    return result;
}

static double triangle_area(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c) {
    geomlib_vec3d ab = geomlib_vec3d_sub(tetra_point(b), tetra_point(a));
    geomlib_vec3d ac = geomlib_vec3d_sub(tetra_point(c), tetra_point(a));
    return 0.5 * geomlib_vec3d_length(geomlib_vec3d_cross(ab, ac));
}

/* In-plane circumcenter; inside/on the closed triangle when *inside is true. */
static geomlib_vec3 triangle_circumcenter(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c,
                                          bool *inside) {
    geomlib_vec3d origin = tetra_point(a);
    geomlib_vec3d u = geomlib_vec3d_sub(tetra_point(b), origin);
    geomlib_vec3d v = geomlib_vec3d_sub(tetra_point(c), origin);
    double u2 = geomlib_vec3d_dot(u, u);
    double v2 = geomlib_vec3d_dot(v, v);
    double uv = geomlib_vec3d_dot(u, v);
    double det = u2 * v2 - uv * uv;
    geomlib_vec3d centroid = geomlib_vec3d_mul(
        geomlib_vec3d_add(geomlib_vec3d_add(tetra_point(a), tetra_point(b)), tetra_point(c)),
        1.0 / 3.0);
    *inside = false;
    if (!(fabs(det) > 1e-30)) return tetra_point_float(centroid);
    double alpha = v2 * (u2 - uv) / (2.0 * det);
    double beta = u2 * (v2 - uv) / (2.0 * det);
    const double tolerance = 1e-9;
    *inside = alpha >= -tolerance && beta >= -tolerance && alpha + beta <= 1.0 + tolerance;
    geomlib_vec3d center = geomlib_vec3d_add(origin,
        geomlib_vec3d_add(geomlib_vec3d_mul(u, alpha), geomlib_vec3d_mul(v, beta)));
    return tetra_point_float(center);
}

static geomlib_vec3 midpoint_of_longest_edge(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c) {
    geomlib_vec3d points[3] = {tetra_point(a), tetra_point(b), tetra_point(c)};
    int longest = 0;
    double best = -1.0;
    const int pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int i = 0; i < 3; i++) {
        geomlib_vec3d edge = geomlib_vec3d_sub(points[pairs[i][0]], points[pairs[i][1]]);
        double length2 = geomlib_vec3d_dot(edge, edge);
        if (length2 > best) {
            best = length2;
            longest = i;
        }
    }
    geomlib_vec3d sum = geomlib_vec3d_add(points[pairs[longest][0]], points[pairs[longest][1]]);
    return tetra_point_float(geomlib_vec3d_mul(sum, 0.5));
}

/* Face-triangle adjacency for the reference surface (up to three neighbors). */
static uint32_t *surface_face_neighbors(const geomlib_trimesh *surface) {
    if (!surface->triangle_count || surface->triangle_count > SIZE_MAX / (3 * sizeof(uint32_t)))
        return NULL;
    uint32_t *neighbors = malloc(surface->triangle_count * 3 * sizeof(*neighbors));
    uint32_t *counts = calloc(surface->triangle_count, sizeof(*counts));
    if (!neighbors || !counts) {
        free(neighbors);
        free(counts);
        return NULL;
    }
    for (size_t i = 0; i < surface->triangle_count * 3; i++) neighbors[i] = UINT32_MAX;
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(surface);
    if (edges.edge_count && !edges.edges) {
        free(neighbors);
        free(counts);
        return NULL;
    }
    for (size_t i = 0; i < edges.edge_count; i++) {
        uint32_t left = edges.edges[i].face0, right = edges.edges[i].face1;
        if (left >= surface->triangle_count || right >= surface->triangle_count) continue;
        if (counts[left] < 3) neighbors[left * 3 + counts[left]++] = right;
        if (counts[right] < 3) neighbors[right * 3 + counts[right]++] = left;
    }
    geomlib_edge_mesh_free(&edges);
    free(counts);
    return neighbors;
}

/* Every boundary face must be contained in one input triangle and the boundary
 * area must match. This accepts a subdivided but geometrically identical shell. */
static int verify_surface_refinement(const geomlib_trimesh *surface, const geomlib_bvh *bvh,
                                     const geomlib_tetramesh *mesh, double tolerance) {
    if (!mesh->boundary_face_count || !surface->triangle_count || !bvh)
        return GEOMLIB_MESH_ERROR_UNSUPPORTED;
    double surface_area = 0.0;
    for (size_t i = 0; i < surface->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(surface, i, &a, &b, &c);
        surface_area += triangle_area(a, b, c);
    }
    if (!(surface_area > 0.0)) return GEOMLIB_MESH_ERROR_UNSUPPORTED;

    uint32_t *neighbors = surface_face_neighbors(surface);
    if (!neighbors) return GEOMLIB_MESH_ERROR_MEMORY;
    int result = GEOMLIB_MESH_ERROR_UNSUPPORTED;
    double boundary_area = 0.0;
    for (size_t face = 0; face < mesh->boundary_face_count; face++) {
        const uint32_t *indices = &mesh->boundary_faces[face * 3];
        geomlib_vec3 points[3] = {mesh->vertices[indices[0]], mesh->vertices[indices[1]],
                                  mesh->vertices[indices[2]]};
        geomlib_vec3d centroid = geomlib_vec3d_mul(geomlib_vec3d_add(
            geomlib_vec3d_add(tetra_point(points[0]), tetra_point(points[1])),
            tetra_point(points[2])), 1.0 / 3.0);
        geomlib_vec3 closest;
        float distance;
        uint32_t primitive;
        if (!geomlib_bvh_closest_point(bvh, tetra_point_float(centroid), &closest, &distance,
                                       &primitive))
            goto finish;
        if (!((double)distance <= tolerance) || primitive >= surface->triangle_count) goto finish;

        bool contained = false;
        for (int slot = -1; slot < 3 && !contained; slot++) {
            uint32_t candidate = slot < 0 ? primitive : neighbors[primitive * 3 + slot];
            if (candidate >= surface->triangle_count) continue;
            geomlib_vec3 a, b, c;
            geomlib_trimesh_get_triangle_vertices(surface, candidate, &a, &b, &c);
            bool inside = true;
            for (int i = 0; i < 3 && inside; i++) {
                geomlib_vec3 nearest = geomlib_closest_point_triangle(points[i], a, b, c);
                inside = geomlib_vec3_length(geomlib_vec3_sub(points[i], nearest)) <= tolerance;
            }
            contained = inside;
        }
        if (!contained) goto finish;
        boundary_area += triangle_area(points[0], points[1], points[2]);
    }
    if (fabs(boundary_area - surface_area) > surface_area * 1e-6 + tolerance * tolerance)
        goto finish;
    result = GEOMLIB_MESH_SUCCESS;

finish:
    free(neighbors);
    return result;
}

static bool sorted_face_contains(const tetra_face_key *keys, size_t count, tetra_face_key key) {
    return bsearch(&key, keys, count, sizeof(*keys), compare_face_keys) != NULL;
}

static int find_missing_facet(const geomlib_trimesh *surface, bool reverse,
                              const geomlib_tetramesh *mesh, size_t *missing) {
    size_t count = mesh->boundary_face_count;
    if (!count) return 0;
    tetra_face_key *keys = malloc(count * sizeof(*keys));
    if (!keys) return -1;
    for (size_t i = 0; i < count; i++)
        keys[i] = oriented_face(&mesh->boundary_faces[i * 3], false);
    qsort(keys, count, sizeof(*keys), compare_face_keys);
    int found = 0;
    for (size_t i = 0; i < surface->triangle_count; i++) {
        tetra_face_key key = oriented_face(&surface->indices[i * 3], reverse);
        if (sorted_face_contains(keys, count, key)) continue;
        *missing = i;
        found = 1;
        break;
    }
    free(keys);
    return found;
}

static bool point_is_existing(const geomlib_tetramesh *mesh, geomlib_vec3 point, double tolerance) {
    double limit = tolerance * tolerance;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        geomlib_vec3d delta = geomlib_vec3d_sub(tetra_point(mesh->vertices[i]), tetra_point(point));
        if (geomlib_vec3d_dot(delta, delta) <= limit) return true;
    }
    return false;
}

/* Insert Steiner points until the boundary is a subdivision of the input surface. */
static int recover_boundary(const geomlib_trimesh *surface, const geomlib_bvh *surface_bvh,
                            geomlib_tetramesh *candidate, bool reverse, double expected_volume6) {
    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(surface);
    double scale = geomlib_vec3d_length(geomlib_vec3d_sub(tetra_point(bounds.max),
                                                          tetra_point(bounds.min)));
    double tolerance = 1e-4 * fmax(scale, 1.0);
    size_t budget = surface->triangle_count < 256 ? surface->triangle_count : 256;
    double expected = fabs(expected_volume6);
    for (size_t attempt = 0; attempt < budget; attempt++) {
        int status = tag_boundary(candidate);
        if (status != GEOMLIB_MESH_SUCCESS)
            return status == GEOMLIB_MESH_ERROR_MEMORY ? status : GEOMLIB_MESH_ERROR_UNSUPPORTED;
        if (verify_surface_refinement(surface, surface_bvh, candidate, tolerance) ==
            GEOMLIB_MESH_SUCCESS) {
            double volume = 0.0;
            for (size_t i = 0; i < candidate->tetrahedron_count; i++)
                volume += cell_volume6(candidate, &candidate->tetrahedra[i]);
            if (fabs(fabs(volume) - expected) > expected * 1e-8)
                return GEOMLIB_MESH_ERROR_UNSUPPORTED;
            return GEOMLIB_MESH_SUCCESS;
        }
        size_t missing;
        int found = find_missing_facet(surface, reverse, candidate, &missing);
        if (found <= 0) return found < 0 ? GEOMLIB_MESH_ERROR_MEMORY : GEOMLIB_MESH_ERROR_UNSUPPORTED;
        const uint32_t *indices = &surface->indices[missing * 3];
        geomlib_vec3 a = surface->vertices[indices[0]];
        geomlib_vec3 b = surface->vertices[indices[1]];
        geomlib_vec3 c = surface->vertices[indices[2]];
        bool inside;
        geomlib_vec3 point = triangle_circumcenter(a, b, c, &inside);
        if (!inside) point = midpoint_of_longest_edge(a, b, c);
        if (point_is_existing(candidate, point, tolerance)) {
            point = tetra_point_float(geomlib_vec3d_mul(geomlib_vec3d_add(
                geomlib_vec3d_add(tetra_point(a), tetra_point(b)), tetra_point(c)), 1.0 / 3.0));
            if (point_is_existing(candidate, point, tolerance))
                return GEOMLIB_MESH_ERROR_UNSUPPORTED;
        }
        status = insert_interior_vertex(candidate, point);
        if (status != GEOMLIB_MESH_SUCCESS) return GEOMLIB_MESH_ERROR_UNSUPPORTED;
    }
    return GEOMLIB_MESH_ERROR_UNSUPPORTED;
}

static int try_conforming_delaunay(const geomlib_trimesh *surface, const geomlib_vec3 *interior,
                                   size_t interior_count, double signed_volume6,
                                   geomlib_tetramesh *output) {
    size_t count = surface->vertex_count + interior_count;
    if (count > SIZE_MAX / sizeof(geomlib_vec3)) return GEOMLIB_MESH_ERROR_MEMORY;
    geomlib_vec3 *points = malloc(count * sizeof(*points));
    if (!points) return GEOMLIB_MESH_ERROR_MEMORY;
    memcpy(points, surface->vertices, surface->vertex_count * sizeof(*points));
    if (interior_count)
        memcpy(points + surface->vertex_count, interior, interior_count * sizeof(*points));
    geomlib_tetramesh candidate = {0};
    int result = geomlib_tetrahedralize_points(points, count, &candidate);
    free(points);
    if (result != GEOMLIB_MESH_SUCCESS) return result == GEOMLIB_MESH_ERROR_MEMORY
        ? result : GEOMLIB_MESH_ERROR_UNSUPPORTED;

    bool *used = interior_count ? calloc(interior_count, sizeof(bool)) : NULL;
    result = GEOMLIB_MESH_ERROR_MEMORY;
    if (interior_count && !used) goto finish;

    size_t kept = 0;
    double volume = 0;
    for (size_t i = 0; i < candidate.tetrahedron_count; i++) {
        const geomlib_tetrahedron cell = candidate.tetrahedra[i];
        geomlib_vec3d vertex[4], center = {0};
        for (int j = 0; j < 4; j++) {
            vertex[j] = tetra_point(candidate.vertices[cell.vertices[j]]);
            center = geomlib_vec3d_add(center, geomlib_vec3d_mul(vertex[j], 0.25));
        }
        float winding = geomlib_compute_winding_number(surface, tetra_point_float(center));
        if (!(winding * (signed_volume6 > 0 ? 1 : -1) > 0.5)) continue;
        volume += oriented_volume6(vertex[0], vertex[1], vertex[2], vertex[3]);
        candidate.tetrahedra[kept++] = cell;
        for (int j = 0; j < 4; j++) {
            if (cell.vertices[j] >= surface->vertex_count) used[cell.vertices[j] - surface->vertex_count] = true;
        }
    }
    candidate.tetrahedron_count = kept;
    result = GEOMLIB_MESH_ERROR_UNSUPPORTED;
    double expected_volume = fabs(signed_volume6);
    if (!kept || fabs(volume - expected_volume) > expected_volume * 1e-8) goto finish;
    for (size_t i = 0; i < interior_count; i++) {
        if (!used[i]) goto finish;
    }

    result = tag_boundary(&candidate);
    if (result != GEOMLIB_MESH_SUCCESS) {
        if (result != GEOMLIB_MESH_ERROR_MEMORY) result = GEOMLIB_MESH_ERROR_UNSUPPORTED;
        goto finish;
    }
    result = verify_surface_boundary(surface, signed_volume6 < 0, &candidate);
    if (result != GEOMLIB_MESH_SUCCESS) {
        geomlib_bvh surface_bvh;
        if (geomlib_bvh_build(&surface_bvh, surface, GEOMLIB_BVH_QUALITY_MEDIUM) ==
            GEOMLIB_MESH_SUCCESS) {
            geomlib_aabb bounds = geomlib_trimesh_compute_aabb(surface);
            double scale = geomlib_vec3d_length(geomlib_vec3d_sub(tetra_point(bounds.max),
                                                                  tetra_point(bounds.min)));
            double tolerance = 1e-4 * fmax(scale, 1.0);
            result = verify_surface_refinement(surface, &surface_bvh, &candidate, tolerance);
            if (result != GEOMLIB_MESH_SUCCESS)
                result = recover_boundary(surface, &surface_bvh, &candidate,
                                          signed_volume6 < 0, signed_volume6);
            geomlib_bvh_free(&surface_bvh);
        } else {
            result = GEOMLIB_MESH_ERROR_UNSUPPORTED;
        }
    }

finish:
    free(used);
    if (result == GEOMLIB_MESH_SUCCESS) {
        geomlib_tetramesh_free(output);
        *output = candidate;
    } else {
        geomlib_tetramesh_free(&candidate);
    }
    return result;
}

int geomlib_tetrahedralize_constrained(const geomlib_trimesh *surface,
                                        const geomlib_vec3 *interior_points, size_t interior_count,
                                        geomlib_tetramesh *tetramesh) {
    if (!surface || !tetramesh || surface->triangle_count < 4 || surface->vertex_count < 4 ||
        surface->vertex_count >= UINT32_MAX || (interior_count && !interior_points) ||
        interior_count > UINT32_MAX - surface->vertex_count - 1)
        return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_validation_result validation = geomlib_validate_mesh(surface, GEOMLIB_VALIDATE_STANDARD);
    bool valid = validation.is_valid;
    geomlib_validation_result_free(&validation);
    if (!valid) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_manifold_analysis topology = geomlib_analyze_manifold(surface);
    bool closed = topology.is_closed && topology.is_manifold;
    geomlib_manifold_analysis_free(&topology);
    if (!closed) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_self_intersection intersections = geomlib_detect_self_intersections(surface);
    bool intersects = intersections.has_intersections;
    geomlib_self_intersection_free(&intersections);
    if (intersects) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_vec3d reference = tetra_point(surface->vertices[surface->indices[0]]);
    geomlib_vec3d moment = {0};
    double total = 0;
    for (size_t i = 0; i < surface->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(surface, i, &a, &b, &c);
        geomlib_vec3d p = geomlib_vec3d_sub(tetra_point(a), reference);
        geomlib_vec3d q = geomlib_vec3d_sub(tetra_point(b), reference);
        geomlib_vec3d r = geomlib_vec3d_sub(tetra_point(c), reference);
        double volume = geomlib_vec3d_dot(p, geomlib_vec3d_cross(q, r));
        total += volume;
        moment = geomlib_vec3d_add(moment, geomlib_vec3d_mul(
            geomlib_vec3d_add(geomlib_vec3d_add(p, q), r), volume));
    }
    if (total == 0) return GEOMLIB_MESH_ERROR_DEGENERATE;

    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(surface);
    double scale = geomlib_vec3d_length(geomlib_vec3d_sub(tetra_point(bounds.max), tetra_point(bounds.min)));
    double tolerance = 32 * DBL_EPSILON * scale;
    double winding = total > 0 ? 1 : -1;
    geomlib_vec3 seed = tetra_point_float(geomlib_vec3d_add(reference,
        geomlib_vec3d_mul(moment, 0.25 / total)));
    for (size_t i = 0; i < interior_count; i++) {
        geomlib_vec3 point = interior_points[i];
        if (!finite_point(point)) return GEOMLIB_MESH_ERROR_INVALID;
        float distance;
        geomlib_find_closest_point_on_mesh(surface, NULL, point, &distance);
        if (!(distance > tolerance) || fabsf(geomlib_compute_winding_number(surface, point)) <= 0.5f)
            return GEOMLIB_MESH_ERROR_INVALID;
        for (size_t j = 0; j < i; j++) {
            geomlib_vec3 previous = interior_points[j];
            if (point.x == previous.x && point.y == previous.y && point.z == previous.z)
                return GEOMLIB_MESH_ERROR_INVALID;
        }
    }
    size_t used_seed = SIZE_MAX;
    for (size_t i = 0; i < interior_count; i++) {
        if (!supports_surface(surface, winding, interior_points[i], tolerance)) continue;
        used_seed = i;
        seed = interior_points[i];
        break;
    }
    if (used_seed == SIZE_MAX && !supports_surface(surface, winding, seed, tolerance))
        return try_conforming_delaunay(surface, interior_points, interior_count, total, tetramesh);

    geomlib_tetramesh result;
    geomlib_tetramesh_init(&result, surface->vertex_count + interior_count + 1, surface->triangle_count);
    int status = GEOMLIB_MESH_ERROR_MEMORY;
    if (!result.vertices || !result.tetrahedra) goto finish;
    memcpy(result.vertices, surface->vertices, surface->vertex_count * sizeof(geomlib_vec3));
    result.vertex_count = surface->vertex_count;
    uint32_t center = geomlib_tetramesh_add_vertex(&result, seed);
    if (center == UINT32_MAX) goto finish;

    for (size_t i = 0; i < surface->triangle_count; i++) {
        const uint32_t *v = &surface->indices[i * 3];
        if (geomlib_tetramesh_add_tetrahedron(&result, center, v[0], v[winding > 0 ? 1 : 2],
                                             v[winding > 0 ? 2 : 1]) == UINT32_MAX) goto finish;
    }
    for (size_t i = 0; i < interior_count; i++) {
        if (i == used_seed) continue;
        status = insert_interior_vertex(&result, interior_points[i]);
        if (status != GEOMLIB_MESH_SUCCESS) goto finish;
    }
    status = tag_boundary(&result);
    if (status != GEOMLIB_MESH_SUCCESS) goto finish;
    status = verify_surface_boundary(surface, winding < 0, &result);
    if (status != GEOMLIB_MESH_SUCCESS) goto finish;
    geomlib_tetramesh_free(tetramesh);
    *tetramesh = result;
    return GEOMLIB_MESH_SUCCESS;

finish:
    geomlib_tetramesh_free(&result);
    return status;
}

int geomlib_tetrahedralize_surface(const geomlib_trimesh *surface, const geomlib_tetra_params *params,
                                    geomlib_tetramesh *tetramesh) {
    if (!surface || !params || !tetramesh) return GEOMLIB_MESH_ERROR_INVALID;
    if (params->method == GEOMLIB_TETRA_OCTREE || params->refine_boundary ||
        params->max_edge_length > 0) return GEOMLIB_MESH_ERROR_UNSUPPORTED;
    if (params->method != GEOMLIB_TETRA_DELAUNAY && params->method != GEOMLIB_TETRA_ADVANCING_FRONT)
        return GEOMLIB_MESH_ERROR_INVALID;
    return geomlib_tetrahedralize_constrained(surface, NULL, 0, tetramesh);
}

/* ============================================================================
 * Quality Metrics
 * ============================================================================ */

geomlib_tetra_quality geomlib_compute_tetra_quality(
    const geomlib_tetramesh *mesh,
    uint32_t tetra_idx
) {
    geomlib_tetra_quality quality;
    memset(&quality, 0, sizeof(quality));

    geomlib_vec3 points[4];
    if (!cell_coordinates(mesh, tetra_idx, points)) return quality;
    geomlib_vec3 v0 = points[0], v1 = points[1], v2 = points[2], v3 = points[3];

    /* Volume */
    quality.volume = geomlib_tetrahedron_volume(v0, v1, v2, v3);

    /* Edge lengths */
    float edges[6];
    edges[0] = geomlib_vec3_length(geomlib_vec3_sub(v1, v0));
    edges[1] = geomlib_vec3_length(geomlib_vec3_sub(v2, v0));
    edges[2] = geomlib_vec3_length(geomlib_vec3_sub(v3, v0));
    edges[3] = geomlib_vec3_length(geomlib_vec3_sub(v2, v1));
    edges[4] = geomlib_vec3_length(geomlib_vec3_sub(v3, v1));
    edges[5] = geomlib_vec3_length(geomlib_vec3_sub(v3, v2));

    float min_edge = edges[0];
    float max_edge = edges[0];
    for (int i = 1; i < 6; i++) {
        if (edges[i] < min_edge) min_edge = edges[i];
        if (edges[i] > max_edge) max_edge = edges[i];
    }

    quality.edge_ratio = (min_edge > 1e-10f) ? (max_edge / min_edge) : FLT_MAX;

    /* Circumradius and inradius for aspect ratio */
    geomlib_vec3 center;
    float circumradius;
    geomlib_tetrahedron_circumsphere(v0, v1, v2, v3, &center, &circumradius);

    /* Inradius = 3 * volume / surface_area */
    float surface_area = 0.0f;
    geomlib_vec3 faces[4][3] = {
        {v0, v1, v2}, {v0, v1, v3}, {v0, v2, v3}, {v1, v2, v3}
    };

    for (int f = 0; f < 4; f++) {
        geomlib_vec3 e1 = geomlib_vec3_sub(faces[f][1], faces[f][0]);
        geomlib_vec3 e2 = geomlib_vec3_sub(faces[f][2], faces[f][0]);
        geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
        surface_area += geomlib_vec3_length(cross) * 0.5f;
    }

    float inradius = (surface_area > 1e-10f) ? (3.0f * quality.volume / surface_area) : 0.0f;
    quality.aspect_ratio = (inradius > 1e-10f) ? (circumradius / inradius) : FLT_MAX;

    /* Shape measure (normalized, 1 = regular tetrahedron) */
    /* For regular tet: aspect_ratio = 3 */
    quality.shape_measure = (quality.aspect_ratio > 0) ?
        (3.0f / quality.aspect_ratio) : 0.0f;
    if (quality.shape_measure > 1.0f) quality.shape_measure = 1.0f;

    /* Sliver detection (very small volume relative to edge lengths) */
    float avg_edge_cubed = powf((max_edge + min_edge) * 0.5f, 3.0f);
    quality.is_sliver = (quality.volume < avg_edge_cubed * 0.001f);

    if (quality.volume > 0) {
        quality.min_dihedral_angle = 180;
        for (int edge = 0; edge < 6; edge++) {
            float angle = geomlib_tetrahedron_dihedral_angle(v0, v1, v2, v3, edge);
            quality.min_dihedral_angle = fminf(quality.min_dihedral_angle, angle);
            quality.max_dihedral_angle = fmaxf(quality.max_dihedral_angle, angle);
        }
    } else {
        quality.is_sliver = true;
    }

    return quality;
}

/* Compute mesh statistics */
void geomlib_tetramesh_compute_stats(
    const geomlib_tetramesh *mesh,
    geomlib_tetramesh_stats *stats
) {
    if (!mesh || !stats) return;

    memset(stats, 0, sizeof(geomlib_tetramesh_stats));

    stats->vertex_count = mesh->vertex_count;
    stats->tetrahedron_count = mesh->tetrahedron_count;
    stats->boundary_face_count = mesh->boundary_face_count;

    if (mesh->tetrahedron_count == 0) return;

    stats->min_volume = FLT_MAX;
    stats->max_volume = -FLT_MAX;
    stats->total_volume = 0.0f;
    stats->min_quality = FLT_MAX;
    stats->avg_quality = 0.0f;
    stats->worst_aspect_ratio = 0.0f;
    stats->sliver_count = 0;

    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        geomlib_tetra_quality q = geomlib_compute_tetra_quality(mesh, i);

        stats->total_volume += q.volume;
        if (q.volume < stats->min_volume) stats->min_volume = q.volume;
        if (q.volume > stats->max_volume) stats->max_volume = q.volume;

        stats->avg_quality += q.shape_measure;
        if (q.shape_measure < stats->min_quality) stats->min_quality = q.shape_measure;
        if (q.aspect_ratio > stats->worst_aspect_ratio) {
            stats->worst_aspect_ratio = q.aspect_ratio;
        }

        if (q.is_sliver) stats->sliver_count++;
    }

    stats->avg_volume = stats->total_volume / mesh->tetrahedron_count;
    stats->avg_quality /= mesh->tetrahedron_count;

    /* Count boundary vertices */
    if (mesh->is_boundary_vertex) {
        for (size_t i = 0; i < mesh->vertex_count; i++) {
            if (mesh->is_boundary_vertex[i]) {
                stats->boundary_vertices++;
            }
        }
    }
}

/* ============================================================================
 * Boundary Operations
 * ============================================================================ */

static size_t count_boundary_faces(const geomlib_fvm_mesh *fvm) {
    size_t count = 0;
    for (size_t i = 0; i < fvm->face_count; i++)
        count += fvm->faces[i].cell_right == UINT32_MAX;
    return count;
}

static int tag_boundary(geomlib_tetramesh *mesh) {
    geomlib_fvm_mesh fvm = {0};
    int result = geomlib_fvm_build(&fvm, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    size_t count = count_boundary_faces(&fvm);
    bool *vertices = mesh->vertex_count ? calloc(mesh->vertex_count, sizeof(bool)) : NULL;
    uint32_t *faces = count ? malloc(count * 3 * sizeof(uint32_t)) : NULL;
    if ((mesh->vertex_count && !vertices) || (count && !faces)) {
        free(vertices);
        free(faces);
        geomlib_fvm_free(&fvm);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    size_t output = 0;
    for (size_t i = 0; i < fvm.face_count; i++) {
        const geomlib_fvm_face *face = &fvm.faces[i];
        if (face->cell_right != UINT32_MAX) continue;
        for (int j = 0; j < 3; j++) {
            uint32_t vertex = face->vertices[j];
            vertices[vertex] = true;
            faces[output++] = vertex;
        }
    }
    free(mesh->boundary_faces);
    free(mesh->is_boundary_vertex);
    mesh->boundary_faces = faces;
    mesh->boundary_face_count = count;
    mesh->is_boundary_vertex = vertices;
    geomlib_fvm_free(&fvm);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_tetramesh_extract_boundary(const geomlib_tetramesh *mesh, geomlib_trimesh *boundary) {
    if (!boundary) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_fvm_mesh fvm = {0};
    int result = geomlib_fvm_build(&fvm, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;
    size_t count = count_boundary_faces(&fvm);
    size_t base = boundary->vertex_count;
    size_t triangles = boundary->triangle_count;
    if (base > UINT32_MAX - mesh->vertex_count ||
        base + mesh->vertex_count > SIZE_MAX / sizeof(geomlib_vec3) ||
        triangles > SIZE_MAX / (3 * sizeof(uint32_t)) - count) {
        result = GEOMLIB_MESH_ERROR_INVALID;
        goto finish;
    }
    result = geomlib_trimesh_reserve_vertices(boundary, base + mesh->vertex_count);
    if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    result = geomlib_trimesh_reserve_triangles(boundary, triangles + count);
    if (result != GEOMLIB_MESH_SUCCESS) goto finish;

    if (mesh->vertex_count) {
        memcpy(boundary->vertices + base, mesh->vertices, mesh->vertex_count * sizeof(geomlib_vec3));
        if (boundary->normals)
            memset(boundary->normals + base, 0, mesh->vertex_count * sizeof(geomlib_vec3));
        if (boundary->texcoords)
            memset(boundary->texcoords + base, 0, mesh->vertex_count * sizeof(geomlib_vec2));
    }
    for (size_t i = 0; i < fvm.face_count; i++) {
        const geomlib_fvm_face *face = &fvm.faces[i];
        if (face->cell_right != UINT32_MAX) continue;
        for (int j = 0; j < 3; j++)
            boundary->indices[triangles * 3 + j] = (uint32_t)base + face->vertices[j];
        triangles++;
    }
    boundary->vertex_count = base + mesh->vertex_count;
    boundary->triangle_count = triangles;

finish:
    geomlib_fvm_free(&fvm);
    return result;
}

static int copy_tetrahedra(const geomlib_tetramesh *source, geomlib_tetramesh *copy) {
    geomlib_tetramesh_init(copy, source->vertex_count, source->tetrahedron_count);
    if ((source->vertex_count && !copy->vertices) ||
        (source->tetrahedron_count && !copy->tetrahedra)) {
        geomlib_tetramesh_free(copy);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    copy->vertex_count = source->vertex_count;
    copy->tetrahedron_count = source->tetrahedron_count;
    if (copy->vertex_count)
        memcpy(copy->vertices, source->vertices, copy->vertex_count * sizeof(*copy->vertices));
    if (copy->tetrahedron_count)
        memcpy(copy->tetrahedra, source->tetrahedra, copy->tetrahedron_count * sizeof(*copy->tetrahedra));
    return GEOMLIB_MESH_SUCCESS;
}

static double cell_volume6(const geomlib_tetramesh *mesh, const geomlib_tetrahedron *cell) {
    const uint32_t *v = cell->vertices;
    return oriented_volume6(tetra_point(mesh->vertices[v[0]]), tetra_point(mesh->vertices[v[1]]),
                            tetra_point(mesh->vertices[v[2]]), tetra_point(mesh->vertices[v[3]]));
}

static bool contains_edge(const geomlib_tetrahedron *cell, uint32_t a, uint32_t b) {
    bool has_a = false, has_b = false;
    for (int i = 0; i < 4; i++) {
        has_a |= cell->vertices[i] == a;
        has_b |= cell->vertices[i] == b;
    }
    return has_a && has_b;
}

static uint64_t edge_key(uint32_t a, uint32_t b) {
    if (a > b) { uint32_t temp = a; a = b; b = temp; }
    return ((uint64_t)a << 32) | b;
}

static int compare_edges(const void *left, const void *right) {
    uint64_t a = *(const uint64_t*)left, b = *(const uint64_t*)right;
    return (a > b) - (a < b);
}

static geomlib_vec3 midpoint(const geomlib_tetramesh *mesh, uint32_t a, uint32_t b) {
    geomlib_vec3d sum = geomlib_vec3d_add(tetra_point(mesh->vertices[a]), tetra_point(mesh->vertices[b]));
    return tetra_point_float(geomlib_vec3d_mul(sum, 0.5));
}

static bool same_point(geomlib_vec3 a, geomlib_vec3 b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static bool longest_refinable_edge(const geomlib_tetramesh *mesh, const uint64_t *locked,
                                    size_t lock_count, double threshold, uint32_t edge[2]) {
    double longest = threshold * threshold;
    bool found = false;
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        const uint32_t *v = mesh->tetrahedra[i].vertices;
        for (int j = 0; j < 4; j++) {
            for (int k = j + 1; k < 4; k++) {
                uint32_t a = v[j], b = v[k];
                geomlib_vec3d delta = geomlib_vec3d_sub(tetra_point(mesh->vertices[a]),
                                                       tetra_point(mesh->vertices[b]));
                double length = geomlib_vec3d_dot(delta, delta);
                if (length <= longest) continue;
                uint64_t key = edge_key(a, b);
                if (lock_count && bsearch(&key, locked, lock_count, sizeof(*locked), compare_edges))
                    continue;
                geomlib_vec3 middle = midpoint(mesh, a, b);
                if (same_point(middle, mesh->vertices[a]) || same_point(middle, mesh->vertices[b])) continue;
                longest = length;
                edge[0] = a;
                edge[1] = b;
                found = true;
            }
        }
    }
    return found;
}

static int split_tetrahedral_edge(geomlib_tetramesh *mesh, uint32_t a, uint32_t b) {
    uint32_t vertex = geomlib_tetramesh_add_vertex(mesh, midpoint(mesh, a, b));
    if (vertex == UINT32_MAX) return GEOMLIB_MESH_ERROR_MEMORY;
    size_t count = mesh->tetrahedron_count;
    for (size_t i = 0; i < count; i++) {
        geomlib_tetrahedron cell = mesh->tetrahedra[i];
        if (!contains_edge(&cell, a, b)) continue;
        geomlib_tetrahedron left = cell, right = cell;
        for (int j = 0; j < 4; j++) {
            if (cell.vertices[j] == a) left.vertices[j] = vertex;
            if (cell.vertices[j] == b) right.vertices[j] = vertex;
        }
        double before = cell_volume6(mesh, &cell);
        double l = cell_volume6(mesh, &left) / before, r = cell_volume6(mesh, &right) / before;
        if (!(l > 0 && r > 0) || fabs(l + r - 1) > 64 * FLT_EPSILON)
            return GEOMLIB_MESH_ERROR_DEGENERATE;
        mesh->tetrahedra[i] = left;
        if (geomlib_tetramesh_add_tetrahedron(mesh, right.vertices[0], right.vertices[1],
                                             right.vertices[2], right.vertices[3]) == UINT32_MAX)
            return GEOMLIB_MESH_ERROR_MEMORY;
    }
    return GEOMLIB_MESH_SUCCESS;
}

static int publish_edit(geomlib_tetramesh *mesh, geomlib_tetramesh *copy, int status, bool edited) {
    if (status == GEOMLIB_MESH_SUCCESS && edited) {
        status = tag_boundary(copy);
        if (status == GEOMLIB_MESH_SUCCESS) {
            geomlib_tetramesh_free(mesh);
            *mesh = *copy;
            return status;
        }
    }
    geomlib_tetramesh_free(copy);
    return status;
}

int geomlib_tetramesh_refine(geomlib_tetramesh *mesh, const geomlib_tetra_params *params) {
    if (!mesh || !params || params->max_iterations < 0 || params->max_insertions < 0 ||
        !isfinite(params->max_edge_length) || params->max_edge_length < 0 ||
        !isfinite(params->min_edge_length) || params->min_edge_length < 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    double target = params->max_edge_length > 0 ? params->max_edge_length : params->element_size;
    if (!(target > 0) || !isfinite(target)) return GEOMLIB_MESH_ERROR_INVALID;
    int budget = params->max_iterations < params->max_insertions
        ? params->max_iterations : params->max_insertions;
    if (!budget || !mesh->tetrahedron_count) return GEOMLIB_MESH_SUCCESS;

    geomlib_fvm_mesh fvm = {0};
    int status = geomlib_fvm_build(&fvm, mesh);
    if (status != GEOMLIB_MESH_SUCCESS) return status;
    size_t lock_count = params->preserve_boundary || !params->refine_boundary
        ? count_boundary_faces(&fvm) * 3 : 0;
    uint64_t *locked = lock_count && lock_count <= SIZE_MAX / sizeof(uint64_t)
        ? malloc(lock_count * sizeof(*locked)) : NULL;
    if (lock_count && !locked) {
        geomlib_fvm_free(&fvm);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    if (lock_count) {
        size_t offset = 0;
        for (size_t i = 0; i < fvm.face_count; i++) {
            const geomlib_fvm_face *face = &fvm.faces[i];
            if (face->cell_right != UINT32_MAX) continue;
            for (int j = 0; j < 3; j++)
                locked[offset++] = edge_key(face->vertices[j], face->vertices[(j + 1) % 3]);
        }
        qsort(locked, lock_count, sizeof(*locked), compare_edges);
    }
    geomlib_fvm_free(&fvm);

    geomlib_tetramesh copy;
    status = copy_tetrahedra(mesh, &copy);
    int splits = 0;
    double threshold = fmax(target, 2.0 * params->min_edge_length);
    while (status == GEOMLIB_MESH_SUCCESS && splits < budget) {
        uint32_t edge[2];
        if (!longest_refinable_edge(&copy, locked, lock_count, threshold, edge)) break;
        status = split_tetrahedral_edge(&copy, edge[0], edge[1]);
        splits++;
    }
    free(locked);
    return publish_edit(mesh, &copy, status, splits != 0);
}

typedef struct flip_quality {
    double worst;
    unsigned slivers;
} flip_quality;

static void record_flip_quality(const geomlib_tetramesh *mesh, const geomlib_tetrahedron *cell,
                                 double volume6, flip_quality *quality) {
    double sum = 0, shortest = DBL_MAX, longest = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            geomlib_vec3d edge = geomlib_vec3d_sub(tetra_point(mesh->vertices[cell->vertices[i]]),
                                                  tetra_point(mesh->vertices[cell->vertices[j]]));
            double length2 = geomlib_vec3d_dot(edge, edge);
            sum += length2;
            shortest = fmin(shortest, length2);
            longest = fmax(longest, length2);
        }
    }
    quality->worst = fmin(quality->worst, fabs(volume6) / (sum * sqrt(sum)));
    double scale = 0.5 * (sqrt(shortest) + sqrt(longest));
    quality->slivers += fabs(volume6) < 0.006 * scale * scale * scale;
}

static uint32_t opposite_corner(const geomlib_tetrahedron *cell, const geomlib_fvm_face *face) {
    for (int i = 0; i < 4; i++) {
        uint32_t vertex = cell->vertices[i];
        if (vertex != face->vertices[0] && vertex != face->vertices[1] && vertex != face->vertices[2])
            return vertex;
    }
    return UINT32_MAX;
}

static int improve_tetrahedral_face(geomlib_tetramesh *mesh, const geomlib_fvm_face *face,
                                    bool slivers_only, bool *edited) {
    if (face->cell_right == UINT32_MAX) return GEOMLIB_MESH_SUCCESS;
    const geomlib_tetrahedron *left = &mesh->tetrahedra[face->cell_left];
    const geomlib_tetrahedron *right = &mesh->tetrahedra[face->cell_right];
    uint32_t a = opposite_corner(left, face), b = opposite_corner(right, face);
    double left_volume = fabs(cell_volume6(mesh, left)), right_volume = fabs(cell_volume6(mesh, right));
    flip_quality old = {DBL_MAX, 0}, replacement = {DBL_MAX, 0};
    record_flip_quality(mesh, left, left_volume, &old);
    record_flip_quality(mesh, right, right_volume, &old);
    if (slivers_only && !old.slivers) return GEOMLIB_MESH_SUCCESS;

    geomlib_tetrahedron cells[3];
    double total = 0;
    for (int i = 0; i < 3; i++) {
        cells[i] = (geomlib_tetrahedron){{a, b, face->vertices[i], face->vertices[(i + 1) % 3]}};
        double volume = cell_volume6(mesh, &cells[i]);
        /* A positive volume on every side places the new edge inside the bipyramid. */
        if (!(volume > 0)) return GEOMLIB_MESH_SUCCESS;
        total += volume;
        record_flip_quality(mesh, &cells[i], volume, &replacement);
    }
    if (fabs(total - left_volume - right_volume) > total * 64 * DBL_EPSILON ||
        replacement.slivers > old.slivers ||
        (slivers_only && replacement.slivers == old.slivers) ||
        !(replacement.worst > old.worst * (1 + 64 * DBL_EPSILON))) return GEOMLIB_MESH_SUCCESS;
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        if (contains_edge(&mesh->tetrahedra[i], a, b)) return GEOMLIB_MESH_SUCCESS;
    }

    if (geomlib_tetramesh_add_tetrahedron(mesh, cells[2].vertices[0], cells[2].vertices[1],
                                         cells[2].vertices[2], cells[2].vertices[3]) == UINT32_MAX)
        return GEOMLIB_MESH_ERROR_MEMORY;
    mesh->tetrahedra[face->cell_left] = cells[0];
    mesh->tetrahedra[face->cell_right] = cells[1];
    *edited = true;
    return GEOMLIB_MESH_SUCCESS;
}

/* Replace three cells sharing an interior edge with two cells sharing a face. */
static int improve_tetrahedral_edge(geomlib_tetramesh *mesh, uint32_t a, uint32_t b,
                                    const uint32_t cell_ids[3], bool slivers_only,
                                    bool *edited) {
    geomlib_tetrahedron old[3] = {mesh->tetrahedra[cell_ids[0]],
                                  mesh->tetrahedra[cell_ids[1]],
                                  mesh->tetrahedra[cell_ids[2]]};
    uint32_t opposite[3];
    for (int i = 0; i < 3; i++) {
        opposite[i] = UINT32_MAX;
        for (int j = 0; j < 4; j++) {
            uint32_t vertex = old[i].vertices[j];
            if (vertex != a && vertex != b) {
                opposite[i] = vertex;
                break;
            }
        }
        if (opposite[i] == UINT32_MAX) return GEOMLIB_MESH_SUCCESS;
    }
    uint32_t c = opposite[0], d = opposite[1], e = opposite[2];
    if (c == d || d == e || c == e) return GEOMLIB_MESH_SUCCESS;

    flip_quality before = {DBL_MAX, 0};
    double old_volume = 0.0;
    for (int i = 0; i < 3; i++) {
        double volume = fabs(cell_volume6(mesh, &old[i]));
        old_volume += volume;
        record_flip_quality(mesh, &old[i], volume, &before);
    }
    if (slivers_only && !before.slivers) return GEOMLIB_MESH_SUCCESS;

    geomlib_vec3d pc = tetra_point(mesh->vertices[c]);
    geomlib_vec3d pd = tetra_point(mesh->vertices[d]);
    geomlib_vec3d pe = tetra_point(mesh->vertices[e]);
    geomlib_tetrahedron cells[2];
    flip_quality after = {DBL_MAX, 0};
    double new_volume = 0.0;
    for (int i = 0; i < 2; i++) {
        uint32_t apex = i == 0 ? a : b;
        double volume6 = oriented_volume6(pc, pd, pe, tetra_point(mesh->vertices[apex]));
        if (!(fabs(volume6) > 0)) return GEOMLIB_MESH_SUCCESS;
        if (volume6 > 0)
            cells[i] = (geomlib_tetrahedron){{c, d, e, apex}};
        else
            cells[i] = (geomlib_tetrahedron){{c, e, d, apex}};
        double volume = fabs(volume6);
        new_volume += volume;
        record_flip_quality(mesh, &cells[i], volume, &after);
    }
    if (fabs(new_volume - old_volume) > old_volume * 64 * DBL_EPSILON ||
        after.slivers > before.slivers || (slivers_only && after.slivers == before.slivers) ||
        !(after.worst > before.worst * (1 + 64 * DBL_EPSILON)))
        return GEOMLIB_MESH_SUCCESS;

    tetra_face_key shared = unordered_face(c, d, e);
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        if (i == cell_ids[0] || i == cell_ids[1] || i == cell_ids[2]) continue;
        const uint32_t *v = mesh->tetrahedra[i].vertices;
        for (int j = 0; j < 4; j++) {
            tetra_face_key other = unordered_face(v[(j + 1) % 4], v[(j + 2) % 4],
                                                  v[(j + 3) % 4]);
            if (compare_face_keys(&shared, &other) == 0) return GEOMLIB_MESH_SUCCESS;
        }
    }

    size_t count = mesh->tetrahedron_count;
    size_t write = 0;
    for (size_t i = 0; i < count; i++) {
        if (i == cell_ids[0] || i == cell_ids[1] || i == cell_ids[2]) continue;
        mesh->tetrahedra[write++] = mesh->tetrahedra[i];
    }
    mesh->tetrahedra[write++] = cells[0];
    mesh->tetrahedra[write++] = cells[1];
    mesh->tetrahedron_count = write;
    *edited = true;
    return GEOMLIB_MESH_SUCCESS;
}

typedef struct edge_cell_record {
    uint64_t key;
    uint32_t cell;
} edge_cell_record;

static int compare_edge_cells(const void *left, const void *right) {
    const edge_cell_record *a = left, *b = right;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (a->cell > b->cell) - (a->cell < b->cell);
}

static int try_edge_flip(geomlib_tetramesh *mesh, bool slivers_only, bool *edited) {
    size_t records = mesh->tetrahedron_count * 6;
    if (mesh->tetrahedron_count > SIZE_MAX / (6 * sizeof(edge_cell_record)))
        return GEOMLIB_MESH_ERROR_MEMORY;
    edge_cell_record *list = malloc(records * sizeof(*list));
    if (!list) return GEOMLIB_MESH_ERROR_MEMORY;
    size_t count = 0;
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        const uint32_t *v = mesh->tetrahedra[i].vertices;
        for (int j = 0; j < 4; j++) {
            for (int k = j + 1; k < 4; k++) {
                list[count++] = (edge_cell_record){edge_key(v[j], v[k]), (uint32_t)i};
            }
        }
    }
    qsort(list, count, sizeof(*list), compare_edge_cells);
    int status = GEOMLIB_MESH_SUCCESS;
    for (size_t i = 0; i < count;) {
        size_t end = i + 1;
        while (end < count && list[end].key == list[i].key) end++;
        if (end - i == 3) {
            uint32_t cells[3] = {list[i].cell, list[i + 1].cell, list[i + 2].cell};
            status = improve_tetrahedral_edge(mesh, (uint32_t)(list[i].key >> 32),
                                              (uint32_t)(list[i].key & UINT32_MAX), cells,
                                              slivers_only, edited);
            if (status != GEOMLIB_MESH_SUCCESS || *edited) break;
        }
        i = end;
    }
    free(list);
    return status;
}

static int improve_tetrahedra(geomlib_tetramesh *mesh, size_t iterations, bool slivers_only) {
    if (!mesh) return GEOMLIB_MESH_ERROR_INVALID;
    if (!iterations || !mesh->tetrahedron_count) return GEOMLIB_MESH_SUCCESS;
    geomlib_fvm_mesh fvm = {0};
    int status = geomlib_fvm_build(&fvm, mesh);
    if (status != GEOMLIB_MESH_SUCCESS) return status;

    geomlib_tetramesh copy;
    status = copy_tetrahedra(mesh, &copy);
    bool edited = false;
    for (size_t i = 0; status == GEOMLIB_MESH_SUCCESS && i < iterations; i++) {
        bool flipped = false;
        status = try_edge_flip(&copy, slivers_only, &flipped);
        for (size_t j = 0; status == GEOMLIB_MESH_SUCCESS && !flipped && j < fvm.face_count; j++)
            status = improve_tetrahedral_face(&copy, &fvm.faces[j], slivers_only, &flipped);
        if (!flipped) break;
        edited = true;
        geomlib_fvm_free(&fvm);
        if (i + 1 < iterations) status = geomlib_fvm_build(&fvm, &copy);
    }
    geomlib_fvm_free(&fvm);
    return publish_edit(mesh, &copy, status, edited);
}

int geomlib_tetramesh_optimize(geomlib_tetramesh *mesh, int iterations) {
    if (iterations < 0) return GEOMLIB_MESH_ERROR_INVALID;
    return improve_tetrahedra(mesh, (size_t)iterations, false);
}

int geomlib_tetramesh_remove_slivers(geomlib_tetramesh *mesh) {
    return improve_tetrahedra(mesh, mesh ? mesh->tetrahedron_count : 0, true);
}

static int compare_vertex_ids(const void *left, const void *right) {
    uint32_t a = *(const uint32_t*)left, b = *(const uint32_t*)right;
    return (a > b) - (a < b);
}

static bool local_quality(const geomlib_tetramesh *mesh, const uint32_t *cells, size_t count,
                          flip_quality *quality) {
    for (size_t i = 0; i < count; i++) {
        if (cells[i] >= mesh->tetrahedron_count) return false;
        double volume = cell_volume6(mesh, &mesh->tetrahedra[cells[i]]);
        record_flip_quality(mesh, &mesh->tetrahedra[cells[i]], volume, quality);
    }
    return true;
}

int geomlib_tetramesh_smooth(geomlib_tetramesh *mesh, int iterations) {
    if (!mesh || iterations < 0) return GEOMLIB_MESH_ERROR_INVALID;
    if (!iterations || !mesh->tetrahedron_count || !mesh->vertex_count)
        return GEOMLIB_MESH_SUCCESS;
    int status = tag_boundary(mesh);
    if (status != GEOMLIB_MESH_SUCCESS) return status;

    size_t vertices = mesh->vertex_count, cells = mesh->tetrahedron_count;
    if (cells > SIZE_MAX / (12 * sizeof(uint32_t)) ||
        vertices > SIZE_MAX / (4 * sizeof(uint32_t)))
        return GEOMLIB_MESH_ERROR_MEMORY;
    size_t *cell_counts = calloc(vertices, sizeof(*cell_counts));
    size_t *cell_offsets = calloc(vertices + 1, sizeof(*cell_offsets));
    size_t *neighbor_counts = calloc(vertices, sizeof(*neighbor_counts));
    size_t *neighbor_offsets = calloc(vertices + 1, sizeof(*neighbor_offsets));
    uint32_t *cell_list = malloc(cells * 4 * sizeof(*cell_list));
    uint32_t *neighbors = malloc(cells * 12 * sizeof(*neighbors));
    if (!cell_counts || !cell_offsets || !neighbor_counts || !neighbor_offsets ||
        !cell_list || !neighbors) {
        free(cell_counts);
        free(cell_offsets);
        free(neighbor_counts);
        free(neighbor_offsets);
        free(cell_list);
        free(neighbors);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t i = 0; i < cells; i++) {
        for (int j = 0; j < 4; j++) {
            uint32_t vertex = mesh->tetrahedra[i].vertices[j];
            if (vertex >= vertices) continue;
            cell_counts[vertex]++;
            neighbor_counts[vertex] += 3;
        }
    }
    for (size_t i = 0; i < vertices; i++) {
        cell_offsets[i + 1] = cell_offsets[i] + cell_counts[i];
        neighbor_offsets[i + 1] = neighbor_offsets[i] + neighbor_counts[i];
    }
    size_t *cell_cursor = malloc(vertices * sizeof(*cell_cursor));
    size_t *neighbor_cursor = malloc(vertices * sizeof(*neighbor_cursor));
    if (!cell_cursor || !neighbor_cursor) {
        free(cell_counts);
        free(cell_offsets);
        free(neighbor_counts);
        free(neighbor_offsets);
        free(cell_list);
        free(neighbors);
        free(cell_cursor);
        free(neighbor_cursor);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memcpy(cell_cursor, cell_offsets, vertices * sizeof(*cell_cursor));
    memcpy(neighbor_cursor, neighbor_offsets, vertices * sizeof(*neighbor_cursor));
    for (size_t i = 0; i < cells; i++) {
        const uint32_t *v = mesh->tetrahedra[i].vertices;
        for (int j = 0; j < 4; j++) {
            uint32_t vertex = v[j];
            if (vertex >= vertices) continue;
            cell_list[cell_cursor[vertex]++] = (uint32_t)i;
            for (int k = 0; k < 4; k++) {
                if (k == j) continue;
                neighbors[neighbor_cursor[vertex]++] = v[k];
            }
        }
    }
    free(cell_cursor);
    free(neighbor_cursor);

    for (size_t i = 0; i < vertices; i++) {
        size_t count = neighbor_counts[i];
        uint32_t *start = neighbors + neighbor_offsets[i];
        qsort(start, count, sizeof(*start), compare_vertex_ids);
        size_t unique = 0;
        for (size_t j = 0; j < count; j++) {
            if (!unique || start[unique - 1] != start[j]) start[unique++] = start[j];
        }
        neighbor_counts[i] = unique;
    }

    for (int pass = 0; pass < iterations; pass++) {
        bool changed = false;
        for (uint32_t v = 0; v < vertices; v++) {
            if (mesh->is_boundary_vertex && mesh->is_boundary_vertex[v]) continue;
            size_t incident = cell_counts[v], ring = neighbor_counts[v];
            if (!incident || !ring) continue;

            flip_quality before = {DBL_MAX, 0};
            if (!local_quality(mesh, cell_list + cell_offsets[v], incident, &before)) continue;

            geomlib_vec3d sum = {0};
            for (size_t i = 0; i < ring; i++)
                sum = geomlib_vec3d_add(sum,
                    tetra_point(mesh->vertices[neighbors[neighbor_offsets[v] + i]]));
            geomlib_vec3 old = mesh->vertices[v];
            mesh->vertices[v] = tetra_point_float(geomlib_vec3d_mul(sum, 1.0 / (double)ring));

            bool valid = true;
            for (size_t i = 0; i < incident && valid; i++)
                valid = cell_volume6(mesh, &mesh->tetrahedra[cell_list[cell_offsets[v] + i]]) > 0;
            flip_quality after = {DBL_MAX, 0};
            if (valid) local_quality(mesh, cell_list + cell_offsets[v], incident, &after);
            if (valid && after.slivers <= before.slivers &&
                after.worst > before.worst * (1 + 1e-6)) {
                changed = true;
            } else {
                mesh->vertices[v] = old;
            }
        }
        if (!changed) break;
    }

    free(cell_counts);
    free(cell_offsets);
    free(neighbor_counts);
    free(neighbor_offsets);
    free(cell_list);
    free(neighbors);
    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_tetrahedron_dihedral_angle(geomlib_vec3 v0, geomlib_vec3 v1,
                                         geomlib_vec3 v2, geomlib_vec3 v3, int edge_idx) {
    static const int endpoints[6][4] = {{0,1,2,3}, {0,2,1,3}, {0,3,1,2},
                                         {1,2,0,3}, {1,3,0,2}, {2,3,0,1}};
    if (edge_idx < 0 || edge_idx >= 6) return NAN;
    geomlib_vec3d points[4] = {tetra_point(v0), tetra_point(v1), tetra_point(v2), tetra_point(v3)};
    const int *indices = endpoints[edge_idx];
    geomlib_vec3d origin = points[indices[0]];
    geomlib_vec3d edge = geomlib_vec3d_sub(points[indices[1]], origin);
    geomlib_vec3d first = geomlib_vec3d_cross(edge, geomlib_vec3d_sub(points[indices[2]], origin));
    geomlib_vec3d second = geomlib_vec3d_cross(edge, geomlib_vec3d_sub(points[indices[3]], origin));
    double norm = geomlib_vec3d_length(first) * geomlib_vec3d_length(second);
    if (!(norm > 0) || !isfinite(norm)) return NAN;
    double cosine = geomlib_vec3d_dot(first, second) / norm;
    return (float)(acos(fmax(-1, fmin(1, cosine))) * (180 / M_PI));
}

void geomlib_tetramesh_mark_boundary(geomlib_tetramesh *mesh) {
    if (mesh) tag_boundary(mesh);
}

static bool valid_export(const geomlib_tetramesh *mesh, const char *filename) {
    if (!filename || !mesh || (mesh->vertex_count && !mesh->vertices) ||
        (mesh->tetrahedron_count && !mesh->tetrahedra)) return false;
    for (size_t i = 0; i < mesh->vertex_count; i++)
        if (!finite_point(mesh->vertices[i])) return false;
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        geomlib_vec3 points[4];
        if (!cell_coordinates(mesh, i, points)) return false;
    }
    return true;
}

static void write_coordinate(FILE *file, float value) {
    char buffer[64];
    snprintf(buffer, sizeof(buffer), "%.*g", FLT_DECIMAL_DIG, value);
    const char *decimal = localeconv()->decimal_point;
    char *position = decimal[0] && strcmp(decimal, ".") ? strstr(buffer, decimal) : NULL;
    if (position) fprintf(file, "%.*s.%s", (int)(position - buffer), buffer, position + strlen(decimal));
    else fputs(buffer, file);
}

static void write_point(FILE *file, geomlib_vec3 point) {
    write_coordinate(file, point.x);
    fputc(' ', file);
    write_coordinate(file, point.y);
    fputc(' ', file);
    write_coordinate(file, point.z);
    fputc('\n', file);
}

static int finish_export(FILE *file) {
    bool failed = ferror(file);
    if (fclose(file)) failed = true;
    return failed ? GEOMLIB_MESH_ERROR_INVALID : GEOMLIB_MESH_SUCCESS;
}

int geomlib_tetramesh_export_vtk(const geomlib_tetramesh *mesh, const char *filename) {
    if (!mesh || mesh->vertex_count > INT32_MAX || mesh->tetrahedron_count > INT32_MAX / 5 ||
        !valid_export(mesh, filename)) return GEOMLIB_MESH_ERROR_INVALID;
    FILE *file = fopen(filename, "w");
    if (!file) return GEOMLIB_MESH_ERROR_INVALID;

    fputs("# vtk DataFile Version 3.0\ngeomlib tetrahedral mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n", file);
    fprintf(file, "POINTS %zu float\n", mesh->vertex_count);
    for (size_t i = 0; i < mesh->vertex_count; i++) write_point(file, mesh->vertices[i]);
    fprintf(file, "CELLS %zu %zu\n", mesh->tetrahedron_count, 5 * mesh->tetrahedron_count);
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        const uint32_t *v = mesh->tetrahedra[i].vertices;
        fprintf(file, "4 %u %u %u %u\n", (unsigned)v[0], (unsigned)v[1], (unsigned)v[2], (unsigned)v[3]);
    }
    fprintf(file, "CELL_TYPES %zu\n", mesh->tetrahedron_count);
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) fputs("10\n", file);
    return finish_export(file);
}

int geomlib_tetramesh_export_gmsh(const geomlib_tetramesh *mesh, const char *filename) {
    if (!mesh || mesh->vertex_count > INT32_MAX || mesh->tetrahedron_count > INT32_MAX ||
        !valid_export(mesh, filename)) return GEOMLIB_MESH_ERROR_INVALID;
    FILE *file = fopen(filename, "w");
    if (!file) return GEOMLIB_MESH_ERROR_INVALID;

    fputs("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n", file);
    fprintf(file, "%zu\n", mesh->vertex_count);
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        fprintf(file, "%zu ", i + 1);
        write_point(file, mesh->vertices[i]);
    }
    fprintf(file, "$EndNodes\n$Elements\n%zu\n", mesh->tetrahedron_count);
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        const uint32_t *v = mesh->tetrahedra[i].vertices;
        fprintf(file, "%zu 4 0 %u %u %u %u\n", i + 1, (unsigned)(v[0] + 1),
                (unsigned)(v[1] + 1), (unsigned)(v[2] + 1), (unsigned)(v[3] + 1));
    }
    fputs("$EndElements\n", file);
    return finish_export(file);
}
