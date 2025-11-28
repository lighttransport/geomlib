/*
 * geomlib - Hole Filling Implementation
 */

#include "geomlib/geomlib_holefill.h"
#include "geomlib/geomlib_smoothing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * Hole Detection
 * ============================================================================ */

geomlib_holes geomlib_detect_holes(const geomlib_trimesh *mesh) {
    geomlib_holes result;

    /* Extract boundary loops */
    geomlib_boundary_loops loops = geomlib_extract_boundary_loops(mesh);

    /* Convert loops to holes */
    result.hole_count = loops.loop_count;
    result.holes = (geomlib_hole*)malloc(result.hole_count * sizeof(geomlib_hole));

    for (size_t i = 0; i < loops.loop_count; i++) {
        result.holes[i] = geomlib_classify_hole(mesh, &loops.loops[i]);
    }

    geomlib_boundary_loops_free(&loops);

    return result;
}

void geomlib_holes_free(geomlib_holes *holes) {
    for (size_t i = 0; i < holes->hole_count; i++) {
        free(holes->holes[i].boundary_vertices);
    }
    free(holes->holes);
}

geomlib_vec3 geomlib_hole_centroid(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole
) {
    geomlib_vec3 centroid = geomlib_vec3_create(0, 0, 0);
    for (size_t i = 0; i < hole->vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[hole->boundary_vertices[i]];
        centroid = geomlib_vec3_add(centroid, v);
    }
    return geomlib_vec3_mul(centroid, 1.0f / (float)hole->vertex_count);
}

geomlib_vec3 geomlib_hole_normal(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole
) {
    geomlib_vec3 normal = geomlib_vec3_create(0, 0, 0);

    /* Average edge normals */
    for (size_t i = 0; i < hole->vertex_count; i++) {
        size_t next = (i + 1) % hole->vertex_count;
        geomlib_vec3 v0 = mesh->vertices[hole->boundary_vertices[i]];
        geomlib_vec3 v1 = mesh->vertices[hole->boundary_vertices[next]];
        geomlib_vec3 edge = geomlib_vec3_sub(v1, v0);

        /* Cross with centroid direction */
        geomlib_vec3 to_center = geomlib_vec3_sub(hole->centroid, v0);
        geomlib_vec3 n = geomlib_vec3_cross(edge, to_center);

        normal = geomlib_vec3_add(normal, n);
    }

    return geomlib_vec3_normalize(normal);
}

float geomlib_hole_area_estimate(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole
) {
    /* Simple fan triangulation area estimate */
    float area = 0.0f;
    geomlib_vec3 center = hole->centroid;

    for (size_t i = 0; i < hole->vertex_count; i++) {
        size_t next = (i + 1) % hole->vertex_count;
        geomlib_vec3 v0 = mesh->vertices[hole->boundary_vertices[i]];
        geomlib_vec3 v1 = mesh->vertices[hole->boundary_vertices[next]];

        geomlib_vec3 e0 = geomlib_vec3_sub(v0, center);
        geomlib_vec3 e1 = geomlib_vec3_sub(v1, center);
        geomlib_vec3 cross = geomlib_vec3_cross(e0, e1);

        area += geomlib_vec3_length(cross) * 0.5f;
    }

    return area;
}

bool geomlib_is_hole_planar(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    float threshold
) {
    if (hole->vertex_count < 3) return true;

    geomlib_vec3 plane_origin, plane_normal;
    geomlib_hole_best_fit_plane(mesh, hole, &plane_origin, &plane_normal);

    /* Check distance of all vertices to plane */
    float max_dist = 0.0f;
    for (size_t i = 0; i < hole->vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[hole->boundary_vertices[i]];
        geomlib_vec3 to_point = geomlib_vec3_sub(v, plane_origin);
        float dist = fabsf(geomlib_vec3_dot(to_point, plane_normal));
        if (dist > max_dist) max_dist = dist;
    }

    /* Normalize by hole size */
    float normalized_dist = max_dist / sqrtf(hole->area_estimate);
    return normalized_dist < threshold;
}

void geomlib_hole_best_fit_plane(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    geomlib_vec3 *plane_origin,
    geomlib_vec3 *plane_normal
) {
    /* Use hole centroid as origin */
    *plane_origin = hole->centroid;

    /* Compute covariance matrix */
    float cov[3][3] = {{0}};

    for (size_t i = 0; i < hole->vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[hole->boundary_vertices[i]];
        geomlib_vec3 d = geomlib_vec3_sub(v, hole->centroid);

        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }

    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    /* Simple approximation: use smallest eigenvector */
    /* For now, use hole normal as approximation */
    *plane_normal = hole->normal;
}

geomlib_hole geomlib_classify_hole(
    const geomlib_trimesh *mesh,
    const geomlib_boundary_loop *loop
) {
    geomlib_hole hole;

    /* Copy boundary */
    hole.vertex_count = loop->vertex_count;
    hole.boundary_vertices = (uint32_t*)malloc(loop->vertex_count * sizeof(uint32_t));
    memcpy(hole.boundary_vertices, loop->vertices, loop->vertex_count * sizeof(uint32_t));

    hole.perimeter = loop->length;

    /* Compute properties */
    hole.centroid = geomlib_hole_centroid(mesh, &hole);
    hole.area_estimate = geomlib_hole_area_estimate(mesh, &hole);
    hole.normal = geomlib_hole_normal(mesh, &hole);

    hole.is_planar = geomlib_is_hole_planar(mesh, &hole, 0.1f);
    hole.planarity = hole.is_planar ? 1.0f : 0.5f;  /* Simplified */

    return hole;
}

/* ============================================================================
 * Hole Filling Parameters
 * ============================================================================ */

geomlib_fill_params geomlib_fill_params_default(geomlib_fill_method method) {
    geomlib_fill_params params;
    params.method = method;
    params.smooth_iterations = 10;
    params.smooth_lambda = 0.5f;
    params.max_edge_length = 0.0f;  /* Auto */
    params.min_angle_deg = 15.0f;
    params.preserve_boundary_edges = true;
    params.fair_interior = true;
    return params;
}

/* ============================================================================
 * Ear Clipping Triangulation (2D)
 * ============================================================================ */

static bool is_ear(
    const geomlib_vec2 *points,
    const bool *active,
    size_t point_count,
    size_t prev,
    size_t curr,
    size_t next
) {
    geomlib_vec2 p0 = points[prev];
    geomlib_vec2 p1 = points[curr];
    geomlib_vec2 p2 = points[next];

    /* Check if triangle is CCW (convex at curr) */
    float cross = (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
    if (cross <= 0) return false;

    /* Check if any other point is inside triangle */
    for (size_t i = 0; i < point_count; i++) {
        if (!active[i] || i == prev || i == curr || i == next) continue;

        geomlib_vec2 p = points[i];

        /* Barycentric test */
        float d0 = (p1.x - p0.x) * (p.y - p0.y) - (p1.y - p0.y) * (p.x - p0.x);
        float d1 = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x);
        float d2 = (p0.x - p2.x) * (p.y - p2.y) - (p0.y - p2.y) * (p.x - p2.x);

        if (d0 >= 0 && d1 >= 0 && d2 >= 0) {
            return false;  /* Point inside */
        }
    }

    return true;
}

int geomlib_ear_clipping(
    const geomlib_vec2 *points,
    size_t point_count,
    uint32_t *triangles,
    size_t *triangle_count
) {
    if (point_count < 3) {
        *triangle_count = 0;
        return GEOMLIB_MESH_SUCCESS;
    }

    bool *active = (bool*)malloc(point_count * sizeof(bool));
    for (size_t i = 0; i < point_count; i++) {
        active[i] = true;
    }

    size_t remaining = point_count;
    size_t tri_count = 0;
    size_t current = 0;

    while (remaining > 3) {
        bool found_ear = false;

        for (size_t attempt = 0; attempt < remaining; attempt++) {
            if (!active[current]) {
                current = (current + 1) % point_count;
                continue;
            }

            /* Find prev and next active vertices */
            size_t prev = current;
            do {
                prev = (prev + point_count - 1) % point_count;
            } while (!active[prev]);

            size_t next = current;
            do {
                next = (next + 1) % point_count;
            } while (!active[next]);

            if (is_ear(points, active, point_count, prev, current, next)) {
                /* Emit triangle */
                triangles[tri_count * 3 + 0] = (uint32_t)prev;
                triangles[tri_count * 3 + 1] = (uint32_t)current;
                triangles[tri_count * 3 + 2] = (uint32_t)next;
                tri_count++;

                /* Remove ear */
                active[current] = false;
                remaining--;
                found_ear = true;
                current = next;
                break;
            }

            current = (current + 1) % point_count;
        }

        if (!found_ear) {
            /* Degenerate polygon - force a triangle */
            size_t indices[3];
            size_t idx = 0;
            for (size_t i = 0; i < point_count && idx < 3; i++) {
                if (active[i]) {
                    indices[idx++] = i;
                }
            }
            if (idx == 3) {
                triangles[tri_count * 3 + 0] = (uint32_t)indices[0];
                triangles[tri_count * 3 + 1] = (uint32_t)indices[1];
                triangles[tri_count * 3 + 2] = (uint32_t)indices[2];
                tri_count++;
                active[indices[1]] = false;
                remaining--;
            } else {
                break;
            }
        }
    }

    /* Emit final triangle */
    if (remaining == 3) {
        size_t indices[3];
        size_t idx = 0;
        for (size_t i = 0; i < point_count && idx < 3; i++) {
            if (active[i]) {
                indices[idx++] = i;
            }
        }
        if (idx == 3) {
            triangles[tri_count * 3 + 0] = (uint32_t)indices[0];
            triangles[tri_count * 3 + 1] = (uint32_t)indices[1];
            triangles[tri_count * 3 + 2] = (uint32_t)indices[2];
            tri_count++;
        }
    }

    *triangle_count = tri_count;
    free(active);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Tangent Plane Filling
 * ============================================================================ */

int geomlib_fill_tangent_plane(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
) {
    if (hole->vertex_count < 3) {
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Compute plane basis */
    geomlib_vec3 origin = hole->centroid;
    geomlib_vec3 normal = hole->normal;

    /* Create orthonormal basis */
    geomlib_vec3 tangent1, tangent2;
    if (fabsf(normal.x) < 0.9f) {
        tangent1 = geomlib_vec3_normalize(geomlib_vec3_cross(
            geomlib_vec3_create(1, 0, 0), normal));
    } else {
        tangent1 = geomlib_vec3_normalize(geomlib_vec3_cross(
            geomlib_vec3_create(0, 1, 0), normal));
    }
    tangent2 = geomlib_vec3_cross(normal, tangent1);

    /* Project to 2D */
    geomlib_vec2 *points_2d = (geomlib_vec2*)malloc(hole->vertex_count * sizeof(geomlib_vec2));
    for (size_t i = 0; i < hole->vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[hole->boundary_vertices[i]];
        geomlib_vec3 d = geomlib_vec3_sub(v, origin);
        points_2d[i].x = geomlib_vec3_dot(d, tangent1);
        points_2d[i].y = geomlib_vec3_dot(d, tangent2);
    }

    /* Triangulate in 2D */
    size_t max_triangles = hole->vertex_count - 2;
    uint32_t *local_triangles = (uint32_t*)malloc(max_triangles * 3 * sizeof(uint32_t));
    size_t triangle_count;

    int result = geomlib_ear_clipping(points_2d, hole->vertex_count,
                                       local_triangles, &triangle_count);

    if (result != GEOMLIB_MESH_SUCCESS) {
        free(points_2d);
        free(local_triangles);
        return result;
    }

    /* Add triangles to mesh */
    size_t old_triangle_count = mesh->triangle_count;
    size_t new_triangle_count = old_triangle_count + triangle_count;

    /* Reallocate if needed */
    if (new_triangle_count > mesh->triangle_capacity) {
        mesh->triangle_capacity = new_triangle_count * 2;
        mesh->indices = (uint32_t*)realloc(mesh->indices,
            mesh->triangle_capacity * 3 * sizeof(uint32_t));
    }

    /* Copy triangles with vertex remapping */
    for (size_t i = 0; i < triangle_count; i++) {
        mesh->indices[(old_triangle_count + i) * 3 + 0] =
            hole->boundary_vertices[local_triangles[i * 3 + 0]];
        mesh->indices[(old_triangle_count + i) * 3 + 1] =
            hole->boundary_vertices[local_triangles[i * 3 + 1]];
        mesh->indices[(old_triangle_count + i) * 3 + 2] =
            hole->boundary_vertices[local_triangles[i * 3 + 2]];
    }

    mesh->triangle_count = new_triangle_count;

    free(points_2d);
    free(local_triangles);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Advancing Front Triangulation
 * ============================================================================ */

/* Helper: Find best vertex to form triangle with edge */
static uint32_t find_best_front_vertex(
    const geomlib_trimesh *mesh,
    uint32_t v0, uint32_t v1,
    const bool *is_boundary,
    const bool *is_used,
    size_t vertex_count,
    geomlib_vec3 normal_hint,
    float *best_quality
) {
    geomlib_vec3 p0 = mesh->vertices[v0];
    geomlib_vec3 p1 = mesh->vertices[v1];
    geomlib_vec3 edge = geomlib_vec3_sub(p1, p0);
    float edge_len = geomlib_vec3_length(edge);

    uint32_t best_v = UINT32_MAX;
    *best_quality = -FLT_MAX;

    for (uint32_t v = 0; v < vertex_count; v++) {
        if (v == v0 || v == v1 || !is_boundary[v] || is_used[v]) {
            continue;
        }

        geomlib_vec3 p2 = mesh->vertices[v];

        /* Compute triangle quality metrics */
        geomlib_vec3 e1 = geomlib_vec3_sub(p2, p0);
        geomlib_vec3 e2 = geomlib_vec3_sub(p2, p1);

        /* Check triangle normal orientation */
        geomlib_vec3 tri_normal = geomlib_vec3_cross(edge, e1);
        float dot = geomlib_vec3_dot(tri_normal, normal_hint);
        if (dot < 0.0f) {
            continue;  /* Wrong orientation */
        }

        /* Compute quality: prefer equilateral, penalize long edges */
        float len1 = geomlib_vec3_length(e1);
        float len2 = geomlib_vec3_length(e2);

        /* Aspect ratio quality */
        float min_len = fminf(edge_len, fminf(len1, len2));
        float max_len = fmaxf(edge_len, fmaxf(len1, len2));
        float aspect = (max_len > 1e-6f) ? (min_len / max_len) : 0.0f;

        /* Angle quality (prefer ~60 degrees) */
        float cos_angle = geomlib_vec3_dot(e1, e2) / (len1 * len2 + 1e-10f);
        float angle_quality = 1.0f - fabsf(cos_angle + 0.5f);  /* -0.5 = cos(120deg) */

        /* Combined quality */
        float quality = aspect * 0.7f + angle_quality * 0.3f;

        if (quality > *best_quality) {
            *best_quality = quality;
            best_v = v;
        }
    }

    return best_v;
}

int geomlib_fill_advancing_front(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
) {
    if (hole->vertex_count < 3) {
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Build active front (list of boundary edges) */
    size_t front_capacity = hole->vertex_count * 2;
    uint32_t *front_v0 = (uint32_t*)malloc(front_capacity * sizeof(uint32_t));
    uint32_t *front_v1 = (uint32_t*)malloc(front_capacity * sizeof(uint32_t));
    bool *front_active = (bool*)calloc(front_capacity, sizeof(bool));

    /* Initialize front with boundary edges */
    size_t front_count = 0;
    for (size_t i = 0; i < hole->vertex_count; i++) {
        front_v0[front_count] = hole->boundary_vertices[i];
        front_v1[front_count] = hole->boundary_vertices[(i + 1) % hole->vertex_count];
        front_active[front_count] = true;
        front_count++;
    }

    /* Track which boundary vertices are available */
    bool *is_boundary = (bool*)calloc(mesh->vertex_count, sizeof(bool));
    bool *is_used = (bool*)calloc(mesh->vertex_count, sizeof(bool));

    for (size_t i = 0; i < hole->vertex_count; i++) {
        is_boundary[hole->boundary_vertices[i]] = true;
    }

    /* Advancing front loop */
    size_t max_iterations = hole->vertex_count * 3;
    size_t iteration = 0;

    while (front_count > 0 && iteration < max_iterations) {
        iteration++;

        /* Find best active edge to process */
        uint32_t best_edge = UINT32_MAX;
        float best_quality = -FLT_MAX;

        for (uint32_t i = 0; i < front_count; i++) {
            if (!front_active[i]) continue;

            float quality;
            uint32_t candidate = find_best_front_vertex(
                mesh,
                front_v0[i], front_v1[i],
                is_boundary, is_used,
                mesh->vertex_count,
                hole->normal,
                &quality
            );

            if (candidate != UINT32_MAX && quality > best_quality) {
                best_quality = quality;
                best_edge = i;
            }
        }

        if (best_edge == UINT32_MAX) {
            /* No suitable vertex found - close remaining front */
            for (uint32_t i = 0; i < front_count; i++) {
                if (front_active[i]) {
                    uint32_t v0 = front_v0[i];
                    uint32_t v1 = front_v1[i];

                    /* Find matching edge to close */
                    for (uint32_t j = i + 1; j < front_count; j++) {
                        if (front_active[j] && front_v0[j] == v1) {
                            uint32_t v2 = front_v1[j];
                            geomlib_trimesh_add_triangle(mesh, v0, v1, v2);
                            front_active[i] = false;
                            front_active[j] = false;
                            break;
                        }
                    }
                }
            }
            break;
        }

        /* Add triangle with best edge */
        uint32_t v0 = front_v0[best_edge];
        uint32_t v1 = front_v1[best_edge];

        float quality;
        uint32_t v2 = find_best_front_vertex(
            mesh, v0, v1,
            is_boundary, is_used,
            mesh->vertex_count,
            hole->normal,
            &quality
        );

        geomlib_trimesh_add_triangle(mesh, v0, v1, v2);

        /* Update front: mark current edge inactive */
        front_active[best_edge] = false;

        /* Add new edges v0->v2 and v2->v1 if not closing hole */
        bool found_v0v2 = false;
        bool found_v2v1 = false;

        for (uint32_t i = 0; i < front_count; i++) {
            if (!front_active[i]) continue;
            if (front_v0[i] == v0 && front_v1[i] == v2) found_v0v2 = true;
            if (front_v0[i] == v2 && front_v1[i] == v1) found_v2v1 = true;
        }

        if (!found_v0v2 && front_count < front_capacity) {
            front_v0[front_count] = v0;
            front_v1[front_count] = v2;
            front_active[front_count] = true;
            front_count++;
        }

        if (!found_v2v1 && front_count < front_capacity) {
            front_v0[front_count] = v2;
            front_v1[front_count] = v1;
            front_active[front_count] = true;
            front_count++;
        }

        is_used[v2] = true;

        /* Count remaining active edges */
        size_t active_count = 0;
        for (uint32_t i = 0; i < front_count; i++) {
            if (front_active[i]) active_count++;
        }

        if (active_count == 0) break;
    }

    /* Cleanup */
    free(front_v0);
    free(front_v1);
    free(front_active);
    free(is_boundary);
    free(is_used);

    /* Apply smoothing if requested */
    if (params->fair_interior && params->smooth_iterations > 0) {
        geomlib_smooth_hole_interior(mesh, hole, hole->boundary_vertices,
                                       hole->vertex_count, params->smooth_iterations);
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Minimal Surface Filling
 * ============================================================================ */

int geomlib_fill_minimal_surface(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
) {
    /* Fill with tangent plane first */
    int result = geomlib_fill_tangent_plane(mesh, hole, params);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return result;
    }

    /* Apply smoothing if requested */
    if (params->fair_interior && params->smooth_iterations > 0) {
        geomlib_smooth_hole_interior(mesh, hole, hole->boundary_vertices,
                                       hole->vertex_count, params->smooth_iterations);
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Hole Filling Main Functions
 * ============================================================================ */

int geomlib_fill_hole(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
) {
    switch (params->method) {
        case GEOMLIB_FILL_TANGENT_PLANE:
            return geomlib_fill_tangent_plane(mesh, hole, params);

        case GEOMLIB_FILL_ADVANCING_FRONT:
            return geomlib_fill_advancing_front(mesh, hole, params);

        case GEOMLIB_FILL_MINIMAL_SURFACE:
            return geomlib_fill_minimal_surface(mesh, hole, params);

        case GEOMLIB_FILL_NATURAL:
            /* TODO: Implement natural fill */
            return geomlib_fill_minimal_surface(mesh, hole, params);

        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }
}

int geomlib_fill_all_holes(
    geomlib_trimesh *mesh,
    const geomlib_fill_params *params
) {
    /* Detect holes */
    geomlib_holes holes = geomlib_detect_holes(mesh);

    /* Fill each hole */
    for (size_t i = 0; i < holes.hole_count; i++) {
        int result = geomlib_fill_hole(mesh, &holes.holes[i], params);
        if (result != GEOMLIB_MESH_SUCCESS) {
            geomlib_holes_free(&holes);
            return result;
        }
    }

    geomlib_holes_free(&holes);
    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Smoothing
 * ============================================================================ */

int geomlib_smooth_hole_interior(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const uint32_t *interior_vertices,
    size_t interior_count,
    uint32_t iterations
) {
    /* Use existing Laplacian smoothing */
    return geomlib_smooth_laplacian(mesh, 0.5f, iterations);
}

int geomlib_fair_filled_region(
    geomlib_trimesh *mesh,
    const uint32_t *region_vertices,
    size_t vertex_count,
    uint32_t iterations,
    float lambda
) {
    return geomlib_smooth_laplacian(mesh, lambda, iterations);
}
