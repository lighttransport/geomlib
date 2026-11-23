/*
 * geomlib - Collision Detection and Mesh Reduction Implementation
 */

#include "geomlib/geomlib_collision.h"
#include "geomlib/geomlib_decimation.h"
#include "geomlib/geomlib_solver.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static void compute_covariance_matrix(const geomlib_vec3 *points, size_t count,
                                       geomlib_vec3 centroid, float cov[3][3]) {
    /* Initialize covariance matrix */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cov[i][j] = 0.0f;
        }
    }

    /* Accumulate covariance */
    for (size_t i = 0; i < count; i++) {
        float dx = points[i].x - centroid.x;
        float dy = points[i].y - centroid.y;
        float dz = points[i].z - centroid.z;

        cov[0][0] += dx * dx;
        cov[0][1] += dx * dy;
        cov[0][2] += dx * dz;
        cov[1][1] += dy * dy;
        cov[1][2] += dy * dz;
        cov[2][2] += dz * dz;
    }

    /* Symmetry */
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    /* Normalize */
    float scale = 1.0f / count;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cov[i][j] *= scale;
        }
    }
}

/* Power iteration for largest eigenvector */
static void power_iteration(float cov[3][3], geomlib_vec3 *eigenvector, float *eigenvalue) {
    geomlib_vec3 v = geomlib_vec3_create(1.0f, 0.0f, 0.0f);

    for (int iter = 0; iter < 20; iter++) {
        geomlib_vec3 v_new;
        v_new.x = cov[0][0] * v.x + cov[0][1] * v.y + cov[0][2] * v.z;
        v_new.y = cov[1][0] * v.x + cov[1][1] * v.y + cov[1][2] * v.z;
        v_new.z = cov[2][0] * v.x + cov[2][1] * v.y + cov[2][2] * v.z;

        float length = geomlib_vec3_length(v_new);
        if (length > 1e-8f) {
            v = geomlib_vec3_normalize(v_new);
        }
    }

    *eigenvector = v;
    *eigenvalue = cov[0][0] * v.x * v.x + cov[1][1] * v.y * v.y + cov[2][2] * v.z * v.z +
                  2.0f * (cov[0][1] * v.x * v.y + cov[0][2] * v.x * v.z + cov[1][2] * v.y * v.z);
}

/* ============================================================================
 * OBB Computation
 * ============================================================================ */

int geomlib_compute_obb(const geomlib_trimesh *mesh, geomlib_obb *obb_out) {
    if (!mesh || !obb_out || mesh->vertex_count < 3) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Compute centroid */
    geomlib_vec3 centroid = geomlib_trimesh_centroid(mesh);

    /* Compute covariance matrix */
    float cov[3][3];
    compute_covariance_matrix(mesh->vertices, mesh->vertex_count, centroid, cov);

    /* Find principal axes using power iteration */
    geomlib_vec3 axis0;
    float lambda0;
    power_iteration(cov, &axis0, &lambda0);

    /* Ensure axis0 is normalized */
    axis0 = geomlib_vec3_normalize(axis0);

    /* Deflate matrix to find second eigenvector */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            float vi = (i == 0) ? axis0.x : ((i == 1) ? axis0.y : axis0.z);
            float vj = (j == 0) ? axis0.x : ((j == 1) ? axis0.y : axis0.z);
            cov[i][j] -= lambda0 * vi * vj;
        }
    }

    geomlib_vec3 axis1;
    float lambda1;
    power_iteration(cov, &axis1, &lambda1);

    /* Ensure axis1 is normalized and orthogonal to axis0 */
    axis1 = geomlib_vec3_normalize(axis1);
    float dot = geomlib_vec3_dot(axis0, axis1);
    axis1 = geomlib_vec3_sub(axis1, geomlib_vec3_mul(axis0, dot));
    axis1 = geomlib_vec3_normalize(axis1);

    /* Third axis is cross product */
    geomlib_vec3 axis2 = geomlib_vec3_cross(axis0, axis1);
    float axis2_len = geomlib_vec3_length(axis2);
    if (axis2_len > 1e-6f) {
        axis2 = geomlib_vec3_normalize(axis2);
    } else {
        /* Fallback: use arbitrary orthogonal vector */
        if (fabsf(axis0.x) < 0.9f) {
            axis2 = geomlib_vec3_cross(axis0, geomlib_vec3_create(1.0f, 0.0f, 0.0f));
        } else {
            axis2 = geomlib_vec3_cross(axis0, geomlib_vec3_create(0.0f, 1.0f, 0.0f));
        }
        axis2 = geomlib_vec3_normalize(axis2);
        axis1 = geomlib_vec3_cross(axis2, axis0);
        axis1 = geomlib_vec3_normalize(axis1);
    }

    /* Set OBB center and axes */
    obb_out->center = centroid;
    obb_out->axes[0] = axis0;
    obb_out->axes[1] = axis1;
    obb_out->axes[2] = axis2;

    /* Compute extents by projecting all vertices */
    float min_proj[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
    float max_proj[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[i];
        geomlib_vec3 rel = geomlib_vec3_sub(v, centroid);

        for (int axis = 0; axis < 3; axis++) {
            float proj = geomlib_vec3_dot(rel, obb_out->axes[axis]);
            if (proj < min_proj[axis]) min_proj[axis] = proj;
            if (proj > max_proj[axis]) max_proj[axis] = proj;
        }
    }

    /* Set half-extents */
    obb_out->half_extents.x = (max_proj[0] - min_proj[0]) * 0.5f;
    obb_out->half_extents.y = (max_proj[1] - min_proj[1]) * 0.5f;
    obb_out->half_extents.z = (max_proj[2] - min_proj[2]) * 0.5f;

    /* Adjust center */
    geomlib_vec3 offset = geomlib_vec3_create(
        (max_proj[0] + min_proj[0]) * 0.5f,
        (max_proj[1] + min_proj[1]) * 0.5f,
        (max_proj[2] + min_proj[2]) * 0.5f
    );

    obb_out->center = geomlib_vec3_add(obb_out->center,
        geomlib_vec3_add(geomlib_vec3_add(
            geomlib_vec3_mul(obb_out->axes[0], offset.x),
            geomlib_vec3_mul(obb_out->axes[1], offset.y)),
            geomlib_vec3_mul(obb_out->axes[2], offset.z)));

    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_obb_volume(const geomlib_obb *obb) {
    return 8.0f * obb->half_extents.x * obb->half_extents.y * obb->half_extents.z;
}

bool geomlib_obb_contains_point(const geomlib_obb *obb, geomlib_vec3 point) {
    geomlib_vec3 rel = geomlib_vec3_sub(point, obb->center);

    for (int i = 0; i < 3; i++) {
        float proj = geomlib_vec3_dot(rel, obb->axes[i]);
        float extent = (i == 0) ? obb->half_extents.x : ((i == 1) ? obb->half_extents.y : obb->half_extents.z);
        if (fabsf(proj) > extent) {
            return false;
        }
    }

    return true;
}

/* ============================================================================
 * Convex Hull - Quickhull Algorithm
 * ============================================================================ */

typedef struct hull_triangle {
    uint32_t v[3];
    bool alive;
} hull_triangle;

typedef struct hull_builder {
    const geomlib_vec3 *points;
    size_t point_count;
    hull_triangle *faces;
    size_t face_count;
    size_t face_capacity;
} hull_builder;

static int hull_add_triangle(hull_builder *builder, uint32_t a, uint32_t b, uint32_t c) {
    if (a == b || b == c || a == c) return GEOMLIB_MESH_SUCCESS;
    if (builder->face_count == builder->face_capacity) {
        size_t capacity = builder->face_capacity ? builder->face_capacity * 2 : 16;
        if (capacity > SIZE_MAX / sizeof(hull_triangle)) return GEOMLIB_MESH_ERROR_MEMORY;
        hull_triangle *faces = realloc(builder->faces, capacity * sizeof(*faces));
        if (!faces) return GEOMLIB_MESH_ERROR_MEMORY;
        builder->faces = faces;
        builder->face_capacity = capacity;
    }
    builder->faces[builder->face_count++] = (hull_triangle){{a, b, c}, true};
    return GEOMLIB_MESH_SUCCESS;
}

static geomlib_vec3d hull_face_normal(const geomlib_vec3 *points, const hull_triangle *face) {
    geomlib_vec3d a = {points[face->v[0]].x, points[face->v[0]].y, points[face->v[0]].z};
    geomlib_vec3d b = {points[face->v[1]].x, points[face->v[1]].y, points[face->v[1]].z};
    geomlib_vec3d c = {points[face->v[2]].x, points[face->v[2]].y, points[face->v[2]].z};
    return geomlib_vec3d_cross(geomlib_vec3d_sub(b, a), geomlib_vec3d_sub(c, a));
}

static double hull_face_distance(const geomlib_vec3 *points, const hull_triangle *face,
                                 size_t point) {
    geomlib_vec3d normal = hull_face_normal(points, face);
    double length = geomlib_vec3d_length(normal);
    if (!(length > 1e-30)) return 0.0;
    geomlib_vec3d origin = {points[face->v[0]].x, points[face->v[0]].y, points[face->v[0]].z};
    geomlib_vec3d query = {points[point].x, points[point].y, points[point].z};
    return geomlib_vec3d_dot(normal, geomlib_vec3d_sub(query, origin)) / length;
}

static int compare_hull_edges(const void *left, const void *right) {
    uint64_t a = *(const uint64_t*)left, b = *(const uint64_t*)right;
    return (a > b) - (a < b);
}

static void hull_orient_face(const geomlib_vec3 *positions, hull_triangle *face, geomlib_vec3 centroid) {
    geomlib_vec3d a = {positions[face->v[0]].x, positions[face->v[0]].y, positions[face->v[0]].z};
    geomlib_vec3d b = {positions[face->v[1]].x, positions[face->v[1]].y, positions[face->v[1]].z};
    geomlib_vec3d c = {positions[face->v[2]].x, positions[face->v[2]].y, positions[face->v[2]].z};
    geomlib_vec3d normal = geomlib_vec3d_cross(geomlib_vec3d_sub(b, a), geomlib_vec3d_sub(c, a));
    geomlib_vec3d center = geomlib_vec3d_mul(geomlib_vec3d_add(geomlib_vec3d_add(a, b), c),
                                             1.0 / 3.0);
    geomlib_vec3d outward = geomlib_vec3d_sub(center,
        (geomlib_vec3d){centroid.x, centroid.y, centroid.z});
    if (geomlib_vec3d_dot(normal, outward) < 0.0) {
        uint32_t swap = face->v[1];
        face->v[1] = face->v[2];
        face->v[2] = swap;
    }
}

/* Incremental convex hull: begin with a tetrahedron and expand the horizon for
 * every point outside the current hull. */
int geomlib_compute_convex_hull(const geomlib_trimesh *mesh, geomlib_trimesh *hull_out) {
    if (!mesh || !mesh->vertices || !hull_out || mesh->vertex_count < 4)
        return GEOMLIB_MESH_ERROR_INVALID;
    const geomlib_vec3 *points = mesh->vertices;
    size_t count = mesh->vertex_count;
    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(mesh);
    double scale = geomlib_vec3d_length((geomlib_vec3d){
        bounds.max.x - bounds.min.x, bounds.max.y - bounds.min.y, bounds.max.z - bounds.min.z});
    if (!(scale > 0.0)) return GEOMLIB_MESH_ERROR_DEGENERATE;
    double epsilon = 1e-7 * scale;

    size_t i0 = 0;
    for (size_t i = 1; i < count; i++)
        if (points[i].x < points[i0].x) i0 = i;
    size_t i1 = SIZE_MAX;
    double best = 0.0;
    for (size_t i = 0; i < count; i++) {
        geomlib_vec3d delta = {points[i].x - points[i0].x, points[i].y - points[i0].y,
                               points[i].z - points[i0].z};
        double distance = geomlib_vec3d_dot(delta, delta);
        if (distance > best) { best = distance; i1 = i; }
    }
    if (i1 == SIZE_MAX || best < epsilon * epsilon) return GEOMLIB_MESH_ERROR_DEGENERATE;

    size_t i2 = SIZE_MAX;
    best = 0.0;
    geomlib_vec3d base = {points[i1].x - points[i0].x, points[i1].y - points[i0].y,
                          points[i1].z - points[i0].z};
    for (size_t i = 0; i < count; i++) {
        geomlib_vec3d delta = {points[i].x - points[i0].x, points[i].y - points[i0].y,
                               points[i].z - points[i0].z};
        double area = geomlib_vec3d_length(geomlib_vec3d_cross(base, delta));
        if (area > best) { best = area; i2 = i; }
    }
    if (i2 == SIZE_MAX || best < epsilon * scale) return GEOMLIB_MESH_ERROR_DEGENERATE;

    geomlib_vec3d normal = geomlib_vec3d_cross(base, (geomlib_vec3d){
        points[i2].x - points[i0].x, points[i2].y - points[i0].y, points[i2].z - points[i0].z});
    size_t i3 = SIZE_MAX;
    best = 0.0;
    for (size_t i = 0; i < count; i++) {
        geomlib_vec3d delta = {points[i].x - points[i0].x, points[i].y - points[i0].y,
                               points[i].z - points[i0].z};
        double height = fabs(geomlib_vec3d_dot(normal, delta));
        if (height > best) { best = height; i3 = i; }
    }
    if (i3 == SIZE_MAX || best < epsilon * scale * scale) return GEOMLIB_MESH_ERROR_DEGENERATE;

    hull_builder builder = {points, count, NULL, 0, 0};
    const uint32_t tetra[4] = {(uint32_t)i0, (uint32_t)i1, (uint32_t)i2, (uint32_t)i3};
    const int face_vertices[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for (int f = 0; f < 4; f++) {
        hull_triangle face = {{tetra[face_vertices[f][0]], tetra[face_vertices[f][1]],
                               tetra[face_vertices[f][2]]}, true};
        int opposite = -1;
        for (int j = 0; j < 4; j++) {
            bool used = false;
            for (int k = 0; k < 3; k++) used |= tetra[j] == face.v[k];
            if (!used) { opposite = j; break; }
        }
        if (opposite >= 0 && hull_face_distance(points, &face, tetra[opposite]) > 0.0) {
            uint32_t swap = face.v[1];
            face.v[1] = face.v[2];
            face.v[2] = swap;
        }
        int status = hull_add_triangle(&builder, face.v[0], face.v[1], face.v[2]);
        if (status != GEOMLIB_MESH_SUCCESS) {
            free(builder.faces);
            return status;
        }
    }

    uint64_t *edges = NULL;
    size_t edge_capacity = 0;
    for (size_t i = 0; i < count; i++) {
        if (i == i0 || i == i1 || i == i2 || i == i3) continue;
        size_t base_faces = builder.face_count;
        bool *visible = calloc(base_faces, sizeof(*visible));
        if (!visible) { free(builder.faces); free(edges); return GEOMLIB_MESH_ERROR_MEMORY; }
        size_t visible_count = 0;
        for (size_t f = 0; f < base_faces; f++) {
            if (!builder.faces[f].alive) continue;
            if (hull_face_distance(points, &builder.faces[f], i) > epsilon) {
                visible[f] = true;
                visible_count++;
            }
        }
        if (!visible_count) { free(visible); continue; }

        size_t needed = 0;
        for (size_t f = 0; f < base_faces; f++)
            if (visible[f]) needed += 3;
        if (needed > edge_capacity) {
            uint64_t *grown = realloc(edges, needed * sizeof(*edges));
            if (!grown) {
                free(visible); free(builder.faces); free(edges);
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
            edges = grown;
            edge_capacity = needed;
        }
        size_t edge_count = 0;
        for (size_t f = 0; f < base_faces; f++) {
            if (!visible[f]) continue;
            for (int j = 0; j < 3; j++) {
                uint32_t a = builder.faces[f].v[j];
                uint32_t b = builder.faces[f].v[(j + 1) % 3];
                edges[edge_count++] = ((uint64_t)a << 32) | b;
            }
        }
        qsort(edges, edge_count, sizeof(*edges), compare_hull_edges);

        int status = GEOMLIB_MESH_SUCCESS;
        for (size_t e = 0; e < edge_count && status == GEOMLIB_MESH_SUCCESS; e++) {
            uint64_t key = edges[e];
            uint32_t a = (uint32_t)(key >> 32), b = (uint32_t)(key & UINT32_MAX);
            uint64_t reverse = ((uint64_t)b << 32) | a;
            if (bsearch(&reverse, edges, edge_count, sizeof(*edges), compare_hull_edges))
                continue;  /* interior edge of the visible patch */
            status = hull_add_triangle(&builder, a, b, (uint32_t)i);
        }
        for (size_t f = 0; f < base_faces; f++)
            if (visible[f]) builder.faces[f].alive = false;
        free(visible);
        if (status != GEOMLIB_MESH_SUCCESS) {
            free(builder.faces);
            free(edges);
            return status;
        }
    }
    free(edges);

    double area_epsilon = 1e-12 * scale * scale;
    for (size_t f = 0; f < builder.face_count; f++) {
        if (!builder.faces[f].alive) continue;
        if (geomlib_vec3d_length(hull_face_normal(points, &builder.faces[f])) <= area_epsilon)
            builder.faces[f].alive = false;
    }

    size_t alive_faces = 0;
    for (size_t f = 0; f < builder.face_count; f++)
        if (builder.faces[f].alive) alive_faces++;
    if (!alive_faces) {
        free(builder.faces);
        return GEOMLIB_MESH_ERROR_DEGENERATE;
    }

    geomlib_vec3d centroid = {0};
    size_t used = 0;
    int32_t *remap = malloc(count * sizeof(*remap));
    uint32_t *used_indices = malloc(count * sizeof(*used_indices));
    if (!remap || !used_indices) {
        free(builder.faces); free(remap); free(used_indices);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < count; i++) remap[i] = -1;
    for (size_t f = 0; f < builder.face_count; f++) {
        if (!builder.faces[f].alive) continue;
        for (int j = 0; j < 3; j++) {
            uint32_t v = builder.faces[f].v[j];
            if (remap[v] < 0) {
                remap[v] = (int32_t)used;
                used_indices[used++] = v;
                centroid = geomlib_vec3d_add(centroid,
                    (geomlib_vec3d){points[v].x, points[v].y, points[v].z});
            }
        }
    }
    if (!used) {
        free(builder.faces); free(remap); free(used_indices);
        return GEOMLIB_MESH_ERROR_DEGENERATE;
    }
    centroid = geomlib_vec3d_mul(centroid, 1.0 / (double)used);
    geomlib_vec3 centroid_float = {(float)centroid.x, (float)centroid.y, (float)centroid.z};

    geomlib_trimesh_init(hull_out, used, alive_faces);
    if ((used && !hull_out->vertices) || (alive_faces && !hull_out->indices)) {
        geomlib_trimesh_free(hull_out);
        free(builder.faces); free(remap); free(used_indices);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < used; i++) {
        if (geomlib_trimesh_add_vertex(hull_out, points[used_indices[i]], NULL, NULL) !=
            GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(hull_out);
            free(builder.faces); free(remap); free(used_indices);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    for (size_t f = 0; f < builder.face_count; f++) {
        if (!builder.faces[f].alive) continue;
        hull_triangle face = builder.faces[f];
        hull_orient_face(points, &face, centroid_float);
        if (geomlib_trimesh_add_triangle(hull_out, remap[face.v[0]], remap[face.v[1]],
                                         remap[face.v[2]]) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(hull_out);
            free(builder.faces); free(remap); free(used_indices);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    free(builder.faces);
    free(remap);
    free(used_indices);
    return GEOMLIB_MESH_SUCCESS;
}

bool geomlib_is_convex(const geomlib_trimesh *mesh, float epsilon) {
    if (!mesh || mesh->triangle_count < 4) {
        return false;
    }

    /* Test each face against all vertices */
    for (size_t tri = 0; tri < mesh->triangle_count; tri++) {
        geomlib_vec3 v0, v1, v2;
        geomlib_trimesh_get_triangle_vertices(mesh, tri, &v0, &v1, &v2);

        geomlib_vec3 normal = geomlib_trimesh_triangle_normal(mesh, tri);

        /* All vertices should be on one side of the plane */
        for (size_t i = 0; i < mesh->vertex_count; i++) {
            geomlib_vec3 v = mesh->vertices[i];
            geomlib_vec3 rel = geomlib_vec3_sub(v, v0);
            float dist = geomlib_vec3_dot(rel, normal);

            if (dist > epsilon) {
                return false;  /* Vertex outside */
            }
        }
    }

    return true;
}

/* ============================================================================
 * Mesh Simplification
 * ============================================================================ */

int geomlib_simplify_mesh(geomlib_trimesh *mesh, size_t target_triangles, int preserve_boundary) {
    if (!mesh || target_triangles < 4 || target_triangles >= mesh->triangle_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_triangle_count = target_triangles;
    params.preserve_boundaries = preserve_boundary != 0;
    params.boundary_weight = preserve_boundary ? 1000.0f : 0.0f;
    geomlib_decimation_stats stats;
    return geomlib_decimate_qem(mesh, &params, &stats);
}

int geomlib_simplify_mesh_ratio(geomlib_trimesh *mesh, float reduction_ratio, int preserve_boundary) {
    if (!mesh || reduction_ratio < 0.0f || reduction_ratio > 1.0f) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    size_t target = (size_t)((1.0f - reduction_ratio) * mesh->triangle_count);
    if (target < 4) target = 4;

    return geomlib_simplify_mesh(mesh, target, preserve_boundary);
}

/* ============================================================================
 * Collision Proxy Generation
 * ============================================================================ */

int geomlib_generate_collision_proxy(
    const geomlib_trimesh *mesh,
    geomlib_proxy_strategy strategy,
    size_t target_triangles,
    geomlib_trimesh *proxy_out)
{
    if (!mesh || !proxy_out) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    switch (strategy) {
        case GEOMLIB_PROXY_CONVEX_HULL:
            return geomlib_compute_convex_hull(mesh, proxy_out);

        case GEOMLIB_PROXY_SIMPLIFIED: {
            /* Copy mesh and simplify */
            proxy_out->vertex_count = mesh->vertex_count;
            proxy_out->triangle_count = mesh->triangle_count;
            proxy_out->vertex_capacity = mesh->vertex_count;
            proxy_out->triangle_capacity = mesh->triangle_count;

            proxy_out->vertices = (geomlib_vec3*)malloc(mesh->vertex_count * sizeof(geomlib_vec3));
            proxy_out->indices = (uint32_t*)malloc(mesh->triangle_count * 3 * sizeof(uint32_t));

            if (!proxy_out->vertices || !proxy_out->indices) {
                free(proxy_out->vertices);
                free(proxy_out->indices);
                return GEOMLIB_MESH_ERROR_MEMORY;
            }

            memcpy(proxy_out->vertices, mesh->vertices, mesh->vertex_count * sizeof(geomlib_vec3));
            memcpy(proxy_out->indices, mesh->indices, mesh->triangle_count * 3 * sizeof(uint32_t));

            proxy_out->normals = NULL;
            proxy_out->texcoords = NULL;

            return geomlib_simplify_mesh(proxy_out, target_triangles, 0);
        }

        case GEOMLIB_PROXY_OBB:
        case GEOMLIB_PROXY_AABB:
            return geomlib_generate_box_proxy(mesh, strategy == GEOMLIB_PROXY_OBB, proxy_out);

        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }
}

int geomlib_generate_box_proxy(const geomlib_trimesh *mesh, int use_obb, geomlib_trimesh *box_out) {
    if (!mesh || !box_out) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (use_obb) {
        /* Use OBB */
        geomlib_obb obb;
        int result = geomlib_compute_obb(mesh, &obb);
        if (result != GEOMLIB_MESH_SUCCESS) {
            return result;
        }

        /* Create box mesh */
        float width = obb.half_extents.x * 2.0f;
        float height = obb.half_extents.y * 2.0f;
        float depth = obb.half_extents.z * 2.0f;

        result = geomlib_create_box(box_out, width, height, depth);
        if (result != GEOMLIB_MESH_SUCCESS) {
            return result;
        }

        /* Transform vertices to OBB orientation */
        for (size_t i = 0; i < box_out->vertex_count; i++) {
            geomlib_vec3 v = box_out->vertices[i];
            geomlib_vec3 rotated = geomlib_vec3_add(
                geomlib_vec3_add(
                    geomlib_vec3_mul(obb.axes[0], v.x),
                    geomlib_vec3_mul(obb.axes[1], v.y)),
                geomlib_vec3_mul(obb.axes[2], v.z));
            box_out->vertices[i] = geomlib_vec3_add(rotated, obb.center);
        }
    } else {
        /* Use AABB */
        geomlib_aabb bounds = geomlib_trimesh_compute_aabb(mesh);

        float width = bounds.max.x - bounds.min.x;
        float height = bounds.max.y - bounds.min.y;
        float depth = bounds.max.z - bounds.min.z;

        int result = geomlib_create_box(box_out, width, height, depth);
        if (result != GEOMLIB_MESH_SUCCESS) {
            return result;
        }

        /* Translate to correct position */
        geomlib_vec3 center = geomlib_aabb_center(bounds);
        for (size_t i = 0; i < box_out->vertex_count; i++) {
            box_out->vertices[i] = geomlib_vec3_add(box_out->vertices[i], center);
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Mesh Analysis
 * ============================================================================ */

int geomlib_analyze_mesh_quality(const geomlib_trimesh *mesh, geomlib_mesh_quality *quality_out) {
    if (!mesh || !quality_out || mesh->triangle_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    quality_out->min_edge_length = FLT_MAX;
    quality_out->max_edge_length = 0.0f;
    quality_out->avg_edge_length = 0.0f;
    quality_out->min_triangle_area = FLT_MAX;
    quality_out->max_triangle_area = 0.0f;
    quality_out->degenerate_count = 0;

    size_t edge_count = 0;

    for (size_t tri = 0; tri < mesh->triangle_count; tri++) {
        geomlib_vec3 v0, v1, v2;
        geomlib_trimesh_get_triangle_vertices(mesh, tri, &v0, &v1, &v2);

        /* Edge lengths */
        float e0 = geomlib_vec3_length(geomlib_vec3_sub(v1, v0));
        float e1 = geomlib_vec3_length(geomlib_vec3_sub(v2, v1));
        float e2 = geomlib_vec3_length(geomlib_vec3_sub(v0, v2));

        quality_out->avg_edge_length += e0 + e1 + e2;
        edge_count += 3;

        if (e0 < quality_out->min_edge_length) quality_out->min_edge_length = e0;
        if (e1 < quality_out->min_edge_length) quality_out->min_edge_length = e1;
        if (e2 < quality_out->min_edge_length) quality_out->min_edge_length = e2;

        if (e0 > quality_out->max_edge_length) quality_out->max_edge_length = e0;
        if (e1 > quality_out->max_edge_length) quality_out->max_edge_length = e1;
        if (e2 > quality_out->max_edge_length) quality_out->max_edge_length = e2;

        /* Triangle area */
        float area = geomlib_trimesh_triangle_area(mesh, tri);
        if (area < quality_out->min_triangle_area) quality_out->min_triangle_area = area;
        if (area > quality_out->max_triangle_area) quality_out->max_triangle_area = area;

        /* Degenerate check */
        if (area < 1e-8f) {
            quality_out->degenerate_count++;
        }
    }

    quality_out->avg_edge_length /= edge_count;

    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_estimate_collision_cost(const geomlib_trimesh *mesh) {
    if (!mesh) {
        return 0.0f;
    }

    /* Simple cost model: triangle count * average edge length */
    /* More triangles = more expensive, larger triangles = more expensive */

    geomlib_mesh_quality quality;
    if (geomlib_analyze_mesh_quality(mesh, &quality) != GEOMLIB_MESH_SUCCESS) {
        return (float)mesh->triangle_count;
    }

    return (float)mesh->triangle_count * quality.avg_edge_length;
}
