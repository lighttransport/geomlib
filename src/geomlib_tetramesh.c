/*
 * geomlib - Tetrahedral Mesh Implementation
 * Delaunay tetrahedralization for CFD simulations
 */

#include "geomlib/geomlib_tetramesh.h"
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
    mesh->vertices = malloc(vertex_capacity * sizeof(geomlib_vec3));
    mesh->tetrahedra = malloc(tetrahedron_capacity * sizeof(geomlib_tetrahedron));
    mesh->boundary_faces = NULL;
    mesh->is_boundary_vertex = NULL;

    mesh->vertex_count = 0;
    mesh->tetrahedron_count = 0;
    mesh->boundary_face_count = 0;
    mesh->vertex_capacity = vertex_capacity;
    mesh->tetrahedron_capacity = tetrahedron_capacity;
}

void geomlib_tetramesh_free(geomlib_tetramesh *mesh) {
    if (!mesh) return;

    free(mesh->vertices);
    free(mesh->tetrahedra);
    free(mesh->boundary_faces);
    free(mesh->is_boundary_vertex);

    memset(mesh, 0, sizeof(geomlib_tetramesh));
}

uint32_t geomlib_tetramesh_add_vertex(
    geomlib_tetramesh *mesh,
    geomlib_vec3 position
) {
    if (mesh->vertex_count >= mesh->vertex_capacity) {
        mesh->vertex_capacity *= 2;
        mesh->vertices = realloc(mesh->vertices,
            mesh->vertex_capacity * sizeof(geomlib_vec3));
    }

    mesh->vertices[mesh->vertex_count] = position;
    return mesh->vertex_count++;
}

uint32_t geomlib_tetramesh_add_tetrahedron(
    geomlib_tetramesh *mesh,
    uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3
) {
    if (mesh->tetrahedron_count >= mesh->tetrahedron_capacity) {
        mesh->tetrahedron_capacity *= 2;
        mesh->tetrahedra = realloc(mesh->tetrahedra,
            mesh->tetrahedron_capacity * sizeof(geomlib_tetrahedron));
    }

    geomlib_tetrahedron tet;
    tet.vertices[0] = v0;
    tet.vertices[1] = v1;
    tet.vertices[2] = v2;
    tet.vertices[3] = v3;

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

/* Check if point is inside circumsphere */
static bool point_in_circumsphere(
    geomlib_vec3 point,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
) {
    geomlib_vec3 center;
    float radius;
    geomlib_tetrahedron_circumsphere(v0, v1, v2, v3, &center, &radius);

    if (radius >= FLT_MAX) return false;

    geomlib_vec3 diff = geomlib_vec3_sub(point, center);
    float dist_sq = geomlib_vec3_dot(diff, diff);
    return dist_sq < (radius * radius - 1e-6f);
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

/* Face structure for cavity computation */
typedef struct face {
    uint32_t v[3];
    uint32_t count;  /* How many tets share this face */
} face;

/* Simple face hash for cavity boundary detection */
static uint64_t face_key(uint32_t a, uint32_t b, uint32_t c) {
    /* Sort vertices */
    if (a > b) { uint32_t t = a; a = b; b = t; }
    if (b > c) { uint32_t t = b; b = c; c = t; }
    if (a > b) { uint32_t t = a; a = b; b = t; }

    /* Pack into uint64_t (assuming indices fit in 21 bits each) */
    return ((uint64_t)a << 42) | ((uint64_t)b << 21) | c;
}

/* Bowyer-Watson incremental insertion */
int geomlib_tetrahedralize_points(
    const geomlib_vec3 *points,
    size_t point_count,
    geomlib_tetramesh *tetramesh
) {
    if (!points || point_count < 4 || !tetramesh) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Compute bounding box */
    geomlib_vec3 min_pt = points[0];
    geomlib_vec3 max_pt = points[0];

    for (size_t i = 1; i < point_count; i++) {
        if (points[i].x < min_pt.x) min_pt.x = points[i].x;
        if (points[i].y < min_pt.y) min_pt.y = points[i].y;
        if (points[i].z < min_pt.z) min_pt.z = points[i].z;
        if (points[i].x > max_pt.x) max_pt.x = points[i].x;
        if (points[i].y > max_pt.y) max_pt.y = points[i].y;
        if (points[i].z > max_pt.z) max_pt.z = points[i].z;
    }

    /* Create super-tetrahedron containing all points */
    geomlib_vec3 center = geomlib_vec3_mul(geomlib_vec3_add(min_pt, max_pt), 0.5f);
    geomlib_vec3 size = geomlib_vec3_sub(max_pt, min_pt);
    float max_dim = fmaxf(fmaxf(size.x, size.y), size.z);
    float scale = max_dim * 10.0f;

    geomlib_vec3 super[4] = {
        geomlib_vec3_create(center.x, center.y + scale, center.z),
        geomlib_vec3_create(center.x - scale, center.y - scale, center.z - scale),
        geomlib_vec3_create(center.x + scale, center.y - scale, center.z - scale),
        geomlib_vec3_create(center.x, center.y - scale, center.z + scale)
    };

    /* Initialize mesh with super-tetrahedron */
    uint32_t sv0 = geomlib_tetramesh_add_vertex(tetramesh, super[0]);
    uint32_t sv1 = geomlib_tetramesh_add_vertex(tetramesh, super[1]);
    uint32_t sv2 = geomlib_tetramesh_add_vertex(tetramesh, super[2]);
    uint32_t sv3 = geomlib_tetramesh_add_vertex(tetramesh, super[3]);

    geomlib_tetramesh_add_tetrahedron(tetramesh, sv0, sv1, sv2, sv3);

    /* Insert points incrementally */
    for (size_t i = 0; i < point_count; i++) {
        uint32_t new_vertex = geomlib_tetramesh_add_vertex(tetramesh, points[i]);

        /* Find all tetrahedra whose circumsphere contains the new point (cavity) */
        bool *bad_tet = calloc(tetramesh->tetrahedron_count, sizeof(bool));
        size_t bad_count = 0;

        for (size_t t = 0; t < tetramesh->tetrahedron_count; t++) {
            geomlib_tetrahedron *tet = &tetramesh->tetrahedra[t];
            geomlib_vec3 v0 = tetramesh->vertices[tet->vertices[0]];
            geomlib_vec3 v1 = tetramesh->vertices[tet->vertices[1]];
            geomlib_vec3 v2 = tetramesh->vertices[tet->vertices[2]];
            geomlib_vec3 v3 = tetramesh->vertices[tet->vertices[3]];

            if (point_in_circumsphere(points[i], v0, v1, v2, v3)) {
                bad_tet[t] = true;
                bad_count++;
            }
        }

        if (bad_count == 0) {
            free(bad_tet);
            continue;
        }

        /* Find boundary faces of cavity (faces belonging to exactly 1 bad tet) */
        #define MAX_FACES 10000
        face *cavity_faces = malloc(MAX_FACES * sizeof(face));
        size_t face_count = 0;

        for (size_t t = 0; t < tetramesh->tetrahedron_count; t++) {
            if (!bad_tet[t]) continue;

            geomlib_tetrahedron *tet = &tetramesh->tetrahedra[t];

            /* Four faces of tetrahedron */
            uint32_t faces[4][3] = {
                {tet->vertices[0], tet->vertices[1], tet->vertices[2]},
                {tet->vertices[0], tet->vertices[1], tet->vertices[3]},
                {tet->vertices[0], tet->vertices[2], tet->vertices[3]},
                {tet->vertices[1], tet->vertices[2], tet->vertices[3]}
            };

            for (int f = 0; f < 4; f++) {
                /* Check if this face is already in the list */
                uint64_t key = face_key(faces[f][0], faces[f][1], faces[f][2]);
                bool found = false;

                for (size_t j = 0; j < face_count; j++) {
                    uint64_t existing = face_key(
                        cavity_faces[j].v[0],
                        cavity_faces[j].v[1],
                        cavity_faces[j].v[2]
                    );
                    if (key == existing) {
                        cavity_faces[j].count++;
                        found = true;
                        break;
                    }
                }

                if (!found && face_count < MAX_FACES) {
                    cavity_faces[face_count].v[0] = faces[f][0];
                    cavity_faces[face_count].v[1] = faces[f][1];
                    cavity_faces[face_count].v[2] = faces[f][2];
                    cavity_faces[face_count].count = 1;
                    face_count++;
                }
            }
        }

        /* Remove bad tetrahedra */
        size_t write_idx = 0;
        for (size_t t = 0; t < tetramesh->tetrahedron_count; t++) {
            if (!bad_tet[t]) {
                if (write_idx != t) {
                    tetramesh->tetrahedra[write_idx] = tetramesh->tetrahedra[t];
                }
                write_idx++;
            }
        }
        tetramesh->tetrahedron_count = write_idx;

        /* Create new tetrahedra from cavity boundary faces to new vertex */
        for (size_t f = 0; f < face_count; f++) {
            if (cavity_faces[f].count == 1) {  /* Boundary face */
                geomlib_tetramesh_add_tetrahedron(
                    tetramesh,
                    cavity_faces[f].v[0],
                    cavity_faces[f].v[1],
                    cavity_faces[f].v[2],
                    new_vertex
                );
            }
        }

        free(cavity_faces);
        free(bad_tet);
    }

    /* Remove tetrahedra containing super-tetrahedron vertices */
    size_t write_idx = 0;
    for (size_t t = 0; t < tetramesh->tetrahedron_count; t++) {
        geomlib_tetrahedron *tet = &tetramesh->tetrahedra[t];
        bool has_super = false;

        for (int v = 0; v < 4; v++) {
            if (tet->vertices[v] <= sv3) {
                has_super = true;
                break;
            }
        }

        if (!has_super) {
            /* Remap vertex indices (shift down by 4) */
            for (int v = 0; v < 4; v++) {
                tet->vertices[v] -= 4;
            }
            if (write_idx != t) {
                tetramesh->tetrahedra[write_idx] = *tet;
            }
            write_idx++;
        }
    }
    tetramesh->tetrahedron_count = write_idx;

    /* Remove super-tetrahedron vertices */
    memmove(tetramesh->vertices, tetramesh->vertices + 4,
        (tetramesh->vertex_count - 4) * sizeof(geomlib_vec3));
    tetramesh->vertex_count -= 4;

    return GEOMLIB_MESH_SUCCESS;
}

/* Generate tetrahedral mesh from bounding box */
int geomlib_tetrahedralize_box(
    geomlib_vec3 min_corner,
    geomlib_vec3 max_corner,
    const geomlib_tetra_params *params,
    geomlib_tetramesh *tetramesh
) {
    if (!params || !tetramesh) return GEOMLIB_MESH_ERROR_INVALID;

    /* Calculate number of subdivisions based on element size */
    geomlib_vec3 size = geomlib_vec3_sub(max_corner, min_corner);
    int nx = (int)(size.x / params->element_size) + 1;
    int ny = (int)(size.y / params->element_size) + 1;
    int nz = (int)(size.z / params->element_size) + 1;

    nx = (nx < 2) ? 2 : nx;
    ny = (ny < 2) ? 2 : ny;
    nz = (nz < 2) ? 2 : nz;

    /* Generate grid points */
    size_t point_count = nx * ny * nz;
    geomlib_vec3 *points = malloc(point_count * sizeof(geomlib_vec3));

    size_t idx = 0;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                float u = (float)i / (nx - 1);
                float v = (float)j / (ny - 1);
                float w = (float)k / (nz - 1);

                points[idx++] = geomlib_vec3_create(
                    min_corner.x + u * size.x,
                    min_corner.y + v * size.y,
                    min_corner.z + w * size.z
                );
            }
        }
    }

    /* Tetrahedralize */
    int result = geomlib_tetrahedralize_points(points, point_count, tetramesh);

    free(points);
    return result;
}

/* Tetrahedralize surface mesh (simplified - just create bounding box mesh) */
int geomlib_tetrahedralize_surface(
    const geomlib_trimesh *surface,
    const geomlib_tetra_params *params,
    geomlib_tetramesh *tetramesh
) {
    if (!surface || !params || !tetramesh) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Compute bounding box of surface */
    geomlib_aabb bbox = geomlib_trimesh_compute_aabb(surface);

    /* Expand slightly */
    geomlib_vec3 expand = geomlib_vec3_create(
        params->element_size,
        params->element_size,
        params->element_size
    );
    bbox.min = geomlib_vec3_sub(bbox.min, expand);
    bbox.max = geomlib_vec3_add(bbox.max, expand);

    /* Generate tetrahedral mesh in bounding box */
    return geomlib_tetrahedralize_box(bbox.min, bbox.max, params, tetramesh);
}

/* Simplified constrained Delaunay (uses surface points) */
int geomlib_tetrahedralize_constrained(
    const geomlib_trimesh *boundary,
    const geomlib_vec3 *interior_points,
    size_t interior_count,
    geomlib_tetramesh *tetramesh
) {
    if (!boundary || !tetramesh) return GEOMLIB_MESH_ERROR_INVALID;

    /* Combine boundary vertices with interior points */
    size_t total_points = boundary->vertex_count + interior_count;
    geomlib_vec3 *all_points = malloc(total_points * sizeof(geomlib_vec3));

    memcpy(all_points, boundary->vertices,
        boundary->vertex_count * sizeof(geomlib_vec3));

    if (interior_points && interior_count > 0) {
        memcpy(all_points + boundary->vertex_count, interior_points,
            interior_count * sizeof(geomlib_vec3));
    }

    int result = geomlib_tetrahedralize_points(all_points, total_points, tetramesh);

    free(all_points);
    return result;
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

    if (tetra_idx >= mesh->tetrahedron_count) return quality;

    geomlib_tetrahedron *tet = &mesh->tetrahedra[tetra_idx];
    geomlib_vec3 v0 = mesh->vertices[tet->vertices[0]];
    geomlib_vec3 v1 = mesh->vertices[tet->vertices[1]];
    geomlib_vec3 v2 = mesh->vertices[tet->vertices[2]];
    geomlib_vec3 v3 = mesh->vertices[tet->vertices[3]];

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

    /* Dihedral angles (simplified - just min/max estimates) */
    quality.min_dihedral_angle = 0.0f;
    quality.max_dihedral_angle = 180.0f;

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

int geomlib_tetramesh_extract_boundary(
    const geomlib_tetramesh *tetramesh,
    geomlib_trimesh *boundary
) {
    if (!tetramesh || !boundary) return GEOMLIB_MESH_ERROR_INVALID;

    /* Count face occurrences */
    #define MAX_BOUNDARY_FACES 100000
    typedef struct {
        uint32_t v[3];
        uint32_t count;
    } boundary_face;

    boundary_face *faces = malloc(MAX_BOUNDARY_FACES * sizeof(boundary_face));
    size_t face_count = 0;

    /* Collect all faces */
    for (size_t t = 0; t < tetramesh->tetrahedron_count; t++) {
        geomlib_tetrahedron *tet = &tetramesh->tetrahedra[t];

        uint32_t tet_faces[4][3] = {
            {tet->vertices[0], tet->vertices[1], tet->vertices[2]},
            {tet->vertices[0], tet->vertices[1], tet->vertices[3]},
            {tet->vertices[0], tet->vertices[2], tet->vertices[3]},
            {tet->vertices[1], tet->vertices[2], tet->vertices[3]}
        };

        for (int f = 0; f < 4; f++) {
            uint64_t key = face_key(tet_faces[f][0], tet_faces[f][1], tet_faces[f][2]);
            bool found = false;

            for (size_t i = 0; i < face_count; i++) {
                uint64_t existing = face_key(faces[i].v[0], faces[i].v[1], faces[i].v[2]);
                if (key == existing) {
                    faces[i].count++;
                    found = true;
                    break;
                }
            }

            if (!found && face_count < MAX_BOUNDARY_FACES) {
                faces[face_count].v[0] = tet_faces[f][0];
                faces[face_count].v[1] = tet_faces[f][1];
                faces[face_count].v[2] = tet_faces[f][2];
                faces[face_count].count = 1;
                face_count++;
            }
        }
    }

    /* Copy vertices */
    for (size_t i = 0; i < tetramesh->vertex_count; i++) {
        geomlib_trimesh_add_vertex(boundary, tetramesh->vertices[i], NULL, NULL);
    }

    /* Add boundary faces (count == 1) */
    for (size_t i = 0; i < face_count; i++) {
        if (faces[i].count == 1) {
            geomlib_trimesh_add_triangle(boundary,
                faces[i].v[0], faces[i].v[1], faces[i].v[2]);
        }
    }

    free(faces);
    return GEOMLIB_MESH_SUCCESS;
}

/* Optimization and refinement (simplified stubs) */
int geomlib_tetramesh_optimize(geomlib_tetramesh *mesh, int iterations) {
    /* Stub - would implement Laplacian smoothing for interior vertices */
    (void)mesh;
    (void)iterations;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_tetramesh_remove_slivers(geomlib_tetramesh *mesh) {
    /* Stub - would remove or improve sliver elements */
    (void)mesh;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_tetramesh_refine(
    geomlib_tetramesh *mesh,
    const geomlib_tetra_params *params
) {
    /* Stub - would insert Steiner points to improve quality */
    (void)mesh;
    (void)params;
    return GEOMLIB_MESH_SUCCESS;
}

/* Dihedral angle computation (stub) */
float geomlib_tetrahedron_dihedral_angle(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3,
    int edge_idx
) {
    (void)v0; (void)v1; (void)v2; (void)v3; (void)edge_idx;
    return 60.0f;  /* Placeholder */
}

void geomlib_tetramesh_mark_boundary(geomlib_tetramesh *mesh) {
    /* Stub - would mark boundary vertices */
    (void)mesh;
}

/* Export functions (stubs) */
int geomlib_tetramesh_export_vtk(
    const geomlib_tetramesh *mesh,
    const char *filename
) {
    (void)mesh; (void)filename;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_tetramesh_export_gmsh(
    const geomlib_tetramesh *mesh,
    const char *filename
) {
    (void)mesh; (void)filename;
    return GEOMLIB_MESH_SUCCESS;
}
