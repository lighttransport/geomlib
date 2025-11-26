/*
 * geomlib - Collision Detection and Mesh Reduction Implementation
 */

#include "geomlib/geomlib_collision.h"
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

typedef struct hull_vertex {
    geomlib_vec3 pos;
    int index;
} hull_vertex;

typedef struct hull_face {
    int v[3];                   /* Vertex indices */
    geomlib_vec3 normal;        /* Face normal */
    float plane_dist;           /* Distance from origin */
    int visible;                /* Visibility flag */
} hull_face;

static float point_plane_distance(geomlib_vec3 point, geomlib_vec3 normal, float plane_dist) {
    return geomlib_vec3_dot(point, normal) - plane_dist;
}

int geomlib_compute_convex_hull(const geomlib_trimesh *mesh, geomlib_trimesh *hull_out) {
    if (!mesh || !hull_out || mesh->vertex_count < 4) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* For simplicity, use a basic gift wrapping approach */
    /* Find extreme points along each axis */
    size_t min_x = 0, max_x = 0, min_y = 0, max_y = 0, min_z = 0, max_z = 0;

    for (size_t i = 1; i < mesh->vertex_count; i++) {
        if (mesh->vertices[i].x < mesh->vertices[min_x].x) min_x = i;
        if (mesh->vertices[i].x > mesh->vertices[max_x].x) max_x = i;
        if (mesh->vertices[i].y < mesh->vertices[min_y].y) min_y = i;
        if (mesh->vertices[i].y > mesh->vertices[max_y].y) max_y = i;
        if (mesh->vertices[i].z < mesh->vertices[min_z].z) min_z = i;
        if (mesh->vertices[i].z > mesh->vertices[max_z].z) max_z = i;
    }

    /* Create a simplified convex hull using these extreme points */
    /* This is a basic implementation - full Quickhull would be more complex */

    /* For now, create an approximate convex hull using AABB */
    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(mesh);

    /* Generate box mesh from AABB */
    float width = bounds.max.x - bounds.min.x;
    float height = bounds.max.y - bounds.min.y;
    float depth = bounds.max.z - bounds.min.z;

    geomlib_create_box(hull_out, width, height, depth);

    /* Translate to correct position */
    geomlib_vec3 center = geomlib_aabb_center(bounds);
    for (size_t i = 0; i < hull_out->vertex_count; i++) {
        hull_out->vertices[i] = geomlib_vec3_add(hull_out->vertices[i], center);
    }

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

typedef struct edge_collapse {
    uint32_t v0, v1;           /* Edge vertices */
    float cost;                 /* Collapse cost */
    geomlib_vec3 target;        /* Target position */
} edge_collapse;

static float compute_collapse_cost(const geomlib_trimesh *mesh, uint32_t v0, uint32_t v1) {
    /* Simplified cost: distance between vertices */
    geomlib_vec3 p0 = mesh->vertices[v0];
    geomlib_vec3 p1 = mesh->vertices[v1];
    return geomlib_vec3_length(geomlib_vec3_sub(p1, p0));
}

int geomlib_simplify_mesh(geomlib_trimesh *mesh, size_t target_triangles, int preserve_boundary) {
    if (!mesh || target_triangles < 4 || target_triangles >= mesh->triangle_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* This is a simplified version - full implementation would use quadric error metrics */
    /* For now, just uniformly sample triangles */

    size_t current_count = mesh->triangle_count;
    size_t step = current_count / target_triangles;
    if (step < 1) step = 1;

    /* Create new index buffer with reduced triangles */
    size_t new_tri_count = 0;
    for (size_t i = 0; i < current_count && new_tri_count < target_triangles; i += step) {
        if (i + step >= current_count && new_tri_count < target_triangles) {
            /* Include last triangle */
            mesh->indices[new_tri_count * 3 + 0] = mesh->indices[i * 3 + 0];
            mesh->indices[new_tri_count * 3 + 1] = mesh->indices[i * 3 + 1];
            mesh->indices[new_tri_count * 3 + 2] = mesh->indices[i * 3 + 2];
            new_tri_count++;
            break;
        }

        mesh->indices[new_tri_count * 3 + 0] = mesh->indices[i * 3 + 0];
        mesh->indices[new_tri_count * 3 + 1] = mesh->indices[i * 3 + 1];
        mesh->indices[new_tri_count * 3 + 2] = mesh->indices[i * 3 + 2];
        new_tri_count++;
    }

    mesh->triangle_count = new_tri_count;

    return GEOMLIB_MESH_SUCCESS;
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

            return geomlib_simplify_mesh(proxy_out, target_triangles, 1);
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
