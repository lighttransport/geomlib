/*
 * geomlib - Hole Filling Implementation
 */

#include "geomlib/geomlib_holefill.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * Hole Detection
 * ============================================================================ */

geomlib_holes geomlib_detect_holes(const geomlib_trimesh *mesh) {
    geomlib_holes result = {0};

    /* Extract boundary loops */
    geomlib_boundary_loops loops = geomlib_extract_boundary_loops(mesh);

    /* Convert loops to holes */
    result.hole_count = loops.loop_count;
    result.holes = (geomlib_hole*)malloc(result.hole_count * sizeof(geomlib_hole));

    if (!result.holes) {
        result.hole_count = 0;
    } else {
        for (size_t i = 0; i < loops.loop_count; i++) {
            result.holes[i] = geomlib_classify_hole(mesh, &loops.loops[i]);
            if (!result.holes[i].boundary_vertices) {
                result.hole_count = i;
                geomlib_holes_free(&result);
                result = (geomlib_holes){0};
                break;
            }
        }
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

/* Jacobi eigendecomposition of a symmetric 3x3 matrix. Eigenvector i is the
 * i-th column of vectors; eigenvalues[i] corresponds to it. */
static void jacobi_eigen_3x3(double a[3][3], double eigenvalues[3], double vectors[3][3]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) vectors[i][j] = i == j ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 32; sweep++) {
        int p = 0, q = 1;
        double largest = fabs(a[0][1]);
        if (fabs(a[0][2]) > largest) { largest = fabs(a[0][2]); p = 0; q = 2; }
        if (fabs(a[1][2]) > largest) { largest = fabs(a[1][2]); p = 1; q = 2; }
        if (largest < 1e-15) break;
        double tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        double t = (tau >= 0.0 ? 1.0 : -1.0) /
                   (fabs(tau) + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        for (int k = 0; k < 3; k++) {
            double akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; k++) {
            double apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; k++) {
            double vkp = vectors[k][p], vkq = vectors[k][q];
            vectors[k][p] = c * vkp - s * vkq;
            vectors[k][q] = s * vkp + c * vkq;
        }
    }
    for (int i = 0; i < 3; i++) eigenvalues[i] = a[i][i];
}

void geomlib_hole_best_fit_plane(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    geomlib_vec3 *plane_origin,
    geomlib_vec3 *plane_normal
) {
    /* Use hole centroid as origin */
    *plane_origin = hole->centroid;

    /* Covariance of the rim; the smallest eigenvector is the plane normal. */
    double cov[3][3] = {{0}};
    for (size_t i = 0; i < hole->vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[hole->boundary_vertices[i]];
        geomlib_vec3 d = geomlib_vec3_sub(v, hole->centroid);
        cov[0][0] += (double)d.x * d.x;
        cov[0][1] += (double)d.x * d.y;
        cov[0][2] += (double)d.x * d.z;
        cov[1][1] += (double)d.y * d.y;
        cov[1][2] += (double)d.y * d.z;
        cov[2][2] += (double)d.z * d.z;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    double eigenvalues[3], vectors[3][3];
    jacobi_eigen_3x3(cov, eigenvalues, vectors);
    int smallest = 0;
    for (int i = 1; i < 3; i++)
        if (eigenvalues[i] < eigenvalues[smallest]) smallest = i;
    geomlib_vec3 normal = {(float)vectors[0][smallest], (float)vectors[1][smallest],
                           (float)vectors[2][smallest]};
    float length = geomlib_vec3_length(normal);
    if (!(length > 1e-12f)) {
        *plane_normal = hole->normal;
        return;
    }
    normal = geomlib_vec3_mul(normal, 1.0f / length);
    if (geomlib_vec3_dot(normal, hole->normal) < 0.0f)
        normal = geomlib_vec3_mul(normal, -1.0f);
    *plane_normal = normal;
}

geomlib_hole geomlib_classify_hole(
    const geomlib_trimesh *mesh,
    const geomlib_boundary_loop *loop
) {
    geomlib_hole hole = {0};

    if (!mesh || !loop || !loop->vertices || loop->vertex_count < 3) return hole;

    /* Copy boundary */
    hole.vertex_count = loop->vertex_count;
    hole.boundary_vertices = (uint32_t*)malloc(loop->vertex_count * sizeof(uint32_t));
    if (!hole.boundary_vertices) return (geomlib_hole){0};
    memcpy(hole.boundary_vertices, loop->vertices, loop->vertex_count * sizeof(uint32_t));

    hole.perimeter = loop->length;

    /* Compute properties */
    hole.centroid = geomlib_hole_centroid(mesh, &hole);
    hole.area_estimate = geomlib_hole_area_estimate(mesh, &hole);
    hole.normal = geomlib_hole_normal(mesh, &hole);

    hole.is_planar = geomlib_is_hole_planar(mesh, &hole, 0.1f);
    geomlib_vec3 plane_origin, plane_normal;
    geomlib_hole_best_fit_plane(mesh, &hole, &plane_origin, &plane_normal);
    float max_distance = 0.0f;
    for (size_t i = 0; i < hole.vertex_count; i++) {
        geomlib_vec3 v = mesh->vertices[hole.boundary_vertices[i]];
        float distance = fabsf(geomlib_vec3_dot(geomlib_vec3_sub(v, plane_origin),
                                                plane_normal));
        if (distance > max_distance) max_distance = distance;
    }
    float size = sqrtf(fmaxf(hole.area_estimate, 1e-12f));
    hole.planarity = 1.0f / (1.0f + max_distance / size);

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

static double orient2d(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c) {
    return ((double)b.x - a.x) * ((double)c.y - a.y) -
           ((double)b.y - a.y) * ((double)c.x - a.x);
}

static bool simple_polygon(const geomlib_vec2 *points, size_t count, double *area) {
    *area = 0.0;
    for (size_t i = 0; i < count; i++) {
        geomlib_vec2 a = points[i], b = points[(i + 1) % count];
        if (!isfinite(a.x) || !isfinite(a.y) || (a.x == b.x && a.y == b.y)) return false;
        *area += orient2d(points[0], a, b);
        for (size_t j = i + 2; j < count; j++) {
            if (i == 0 && j == count - 1) continue;
            geomlib_vec2 c = points[j], d = points[(j + 1) % count];
            if (fmaxf(a.x, b.x) < fminf(c.x, d.x) || fmaxf(c.x, d.x) < fminf(a.x, b.x) ||
                fmaxf(a.y, b.y) < fminf(c.y, d.y) || fmaxf(c.y, d.y) < fminf(a.y, b.y)) continue;
            double u = orient2d(a, b, c), v = orient2d(a, b, d);
            double w = orient2d(c, d, a), x = orient2d(c, d, b);
            if (u * v <= 0.0 && w * x <= 0.0) return false;
        }
    }
    return *area != 0.0;
}

static bool is_ear(const geomlib_vec2 *points, const uint32_t *next, uint32_t a,
                   uint32_t b, uint32_t c, double winding) {
    if (winding * orient2d(points[a], points[b], points[c]) <= 0.0) return false;
    for (uint32_t i = next[c]; i != a; i = next[i]) {
        if (winding * orient2d(points[a], points[b], points[i]) >= 0.0 &&
            winding * orient2d(points[b], points[c], points[i]) >= 0.0 &&
            winding * orient2d(points[c], points[a], points[i]) >= 0.0) return false;
    }
    return true;
}

static int triangulate_front(const geomlib_vec2 *points, size_t count, bool best_ear,
                             uint32_t *triangles, size_t *triangle_count) {
    if (!triangle_count) return GEOMLIB_MESH_ERROR_INVALID;
    *triangle_count = 0;
    if (count < 3) return GEOMLIB_MESH_SUCCESS;
    if (!points || !triangles || count >= UINT32_MAX || count > SIZE_MAX / (3 * sizeof(uint32_t)))
        return GEOMLIB_MESH_ERROR_INVALID;
    double area;
    if (!simple_polygon(points, count, &area)) return GEOMLIB_MESH_ERROR_DEGENERATE;
    double winding = area > 0.0 ? 1.0 : -1.0;

    uint32_t *next = (uint32_t*)malloc(count * sizeof(uint32_t));
    if (!next) return GEOMLIB_MESH_ERROR_MEMORY;
    for (size_t i = 0; i < count; i++) next[i] = (uint32_t)((i + 1) % count);

    uint32_t start = 0;
    for (size_t remaining = count; remaining >= 3; remaining--) {
        uint32_t a = start, chosen = UINT32_MAX;
        double best = -1.0;
        do {
            uint32_t b = next[a], c = next[b];
            if (is_ear(points, next, a, b, c, winding)) {
                double dx = (double)points[a].x - points[c].x;
                double dy = (double)points[a].y - points[c].y;
                double quality = fabs(orient2d(points[a], points[b], points[c])) / (dx * dx + dy * dy);
                if (quality > best) {
                    chosen = a;
                    best = quality;
                }
                if (!best_ear) break;
            }
            a = next[a];
        } while (a != start);
        if (chosen == UINT32_MAX) {
            free(next);
            return GEOMLIB_MESH_ERROR_DEGENERATE;
        }
        uint32_t b = next[chosen], c = next[b];
        *triangles++ = chosen;
        *triangles++ = b;
        *triangles++ = c;
        next[chosen] = c;
        start = chosen;
    }
    free(next);
    *triangle_count = count - 2;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_ear_clipping(const geomlib_vec2 *points, size_t point_count, uint32_t *triangles,
                         size_t *triangle_count) {
    return triangulate_front(points, point_count, false, triangles, triangle_count);
}

int geomlib_triangulate_polygon_2d(const geomlib_vec2 *points, size_t point_count,
                                    uint32_t *triangles, size_t *triangle_count) {
    return geomlib_ear_clipping(points, point_count, triangles, triangle_count);
}

/* The cap must traverse each shared rim edge opposite to the original face. */
static int rim_orientation(const geomlib_trimesh *mesh, const geomlib_hole *hole) {
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);
    int orientation = 0;
    for (size_t i = 0; i < hole->vertex_count; i++) {
        uint32_t a = hole->boundary_vertices[i], b = hole->boundary_vertices[(i + 1) % hole->vertex_count];
        if (a >= edges.vertex_count || b >= edges.vertex_count || a == b) {
            orientation = 0;
            break;
        }
        const geomlib_edge *boundary = NULL;
        size_t degree = 0;
        for (uint32_t j = 0; j < edges.vertex_edge_count[a]; j++) {
            const geomlib_edge *edge = &edges.edges[edges.vertex_edges[a][j]];
            if (edge->face1 != UINT32_MAX) continue;
            degree++;
            if (edge->v0 == b || edge->v1 == b) boundary = edge;
        }
        if (!boundary || degree != 2) {
            orientation = 0;
            break;
        }
        const uint32_t *face = &mesh->indices[boundary->face0 * 3];
        int direction = 1;
        for (int k = 0; k < 3; k++) {
            if (face[k] == a && face[(k + 1) % 3] == b) direction = -1;
        }
        if (orientation && orientation != direction) {
            orientation = 0;
            break;
        }
        orientation = direction;
    }
    geomlib_edge_mesh_free(&edges);
    return orientation;
}

/* Split every internal diagonal and join the resulting rings to triangle centers. */
static int refine_patch(geomlib_trimesh *patch) {
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(patch);
    size_t count = patch->triangle_count;
    uint32_t *midpoints = (uint32_t*)malloc(count * 3 * sizeof(uint32_t));
    uint32_t *refined = (uint32_t*)malloc(count * 18 * sizeof(uint32_t));
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (!edges.edge_count || !midpoints || !refined) goto finish;
    memset(midpoints, 0xff, count * 3 * sizeof(uint32_t));
    for (size_t i = 0; i < edges.edge_count; i++) {
        const geomlib_edge *edge = &edges.edges[i];
        if (edge->face1 == UINT32_MAX) continue;
        uint32_t midpoint = (uint32_t)patch->vertex_count;
        geomlib_vec3 point = geomlib_vec3_mul(
            geomlib_vec3_add(patch->vertices[edge->v0], patch->vertices[edge->v1]), 0.5f);
        geomlib_vec2 uv;
        if (patch->texcoords) {
            geomlib_vec2 a = patch->texcoords[edge->v0], b = patch->texcoords[edge->v1];
            uv = (geomlib_vec2){0.5f * (a.x + b.x), 0.5f * (a.y + b.y)};
        }
        if (geomlib_trimesh_add_vertex(patch, point, NULL, patch->texcoords ? &uv : NULL)
            != GEOMLIB_MESH_SUCCESS) goto finish;
        uint32_t faces[] = {edge->face0, edge->face1};
        for (int f = 0; f < 2; f++) {
            uint32_t *indices = &patch->indices[3 * faces[f]];
            for (int j = 0; j < 3; j++) {
                uint32_t a = indices[j], b = indices[(j + 1) % 3];
                if ((a == edge->v0 && b == edge->v1) || (a == edge->v1 && b == edge->v0))
                    midpoints[3 * faces[f] + j] = midpoint;
            }
        }
    }

    size_t face_count = 0;
    for (size_t i = 0; i < count; i++) {
        uint32_t ring[6], length = 0;
        geomlib_vec3 point = {0};
        geomlib_vec2 uv = {0};
        for (int j = 0; j < 3; j++) {
            uint32_t index = patch->indices[3 * i + j];
            point = geomlib_vec3_add(point, geomlib_vec3_mul(patch->vertices[index], 1.0f / 3.0f));
            if (patch->texcoords) {
                uv.x += patch->texcoords[index].x / 3.0f;
                uv.y += patch->texcoords[index].y / 3.0f;
            }
            ring[length++] = index;
            uint32_t midpoint = midpoints[3 * i + j];
            if (midpoint != UINT32_MAX) ring[length++] = midpoint;
        }
        uint32_t center = (uint32_t)patch->vertex_count;
        if (geomlib_trimesh_add_vertex(patch, point, NULL, patch->texcoords ? &uv : NULL)
            != GEOMLIB_MESH_SUCCESS) goto finish;
        for (uint32_t j = 0; j < length; j++) {
            refined[3 * face_count] = ring[j];
            refined[3 * face_count + 1] = ring[(j + 1) % length];
            refined[3 * face_count + 2] = center;
            face_count++;
        }
    }
    free(patch->indices);
    patch->indices = refined;
    patch->triangle_count = face_count;
    patch->triangle_capacity = count * 6;
    refined = NULL;
    result = GEOMLIB_MESH_SUCCESS;

finish:
    geomlib_edge_mesh_free(&edges);
    free(midpoints);
    free(refined);
    return result;
}

static int fill_patch(geomlib_trimesh *mesh, const geomlib_hole *hole,
                      const geomlib_fill_params *params, bool best_ear, bool fair) {
    if (!mesh || !mesh->vertices || !hole || !params || !hole->boundary_vertices ||
        hole->vertex_count < 3 || hole->vertex_count >= UINT32_MAX / 3 ||
        hole->vertex_count > SIZE_MAX / (18 * sizeof(uint32_t)))
        return GEOMLIB_MESH_ERROR_INVALID;
    if (fair && (!isfinite(params->smooth_lambda) || params->smooth_lambda < 0.0f ||
                 params->smooth_lambda > 1.0f)) return GEOMLIB_MESH_ERROR_INVALID;

    size_t count = hole->vertex_count;
    int orientation = rim_orientation(mesh, hole);
    if (!orientation) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_vec3 normal = geomlib_hole_normal(mesh, hole);
    if (!(geomlib_vec3_length(normal) > 0.0f)) return GEOMLIB_MESH_ERROR_DEGENERATE;
    geomlib_vec3 axis = fabsf(normal.x) < 0.9f ? (geomlib_vec3){1, 0, 0} : (geomlib_vec3){0, 1, 0};
    geomlib_vec3 u = geomlib_vec3_normalize(geomlib_vec3_cross(axis, normal));
    geomlib_vec3 v = geomlib_vec3_cross(normal, u);
    geomlib_vec3 origin = mesh->vertices[hole->boundary_vertices[0]];
    geomlib_vec2 *points = (geomlib_vec2*)malloc(count * sizeof(geomlib_vec2));
    if (!points) return GEOMLIB_MESH_ERROR_MEMORY;

    geomlib_trimesh patch = {0};
    uint32_t *interior = NULL;
    int result = geomlib_trimesh_init(&patch, count, count - 2);
    if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    for (size_t i = 0; i < count; i++) {
        uint32_t vertex = hole->boundary_vertices[i];
        geomlib_vec3 point = mesh->vertices[vertex];
        geomlib_vec3 d = geomlib_vec3_sub(point, origin);
        points[i] = (geomlib_vec2){geomlib_vec3_dot(d, u), geomlib_vec3_dot(d, v)};
        result = geomlib_trimesh_add_vertex(&patch, point, NULL,
                                           mesh->texcoords ? &mesh->texcoords[vertex] : NULL);
        if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    }
    result = triangulate_front(points, count, best_ear, patch.indices, &patch.triangle_count);
    if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    if (orientation < 0) {
        for (size_t i = 0; i < patch.triangle_count; i++) {
            uint32_t temp = patch.indices[3 * i + 1];
            patch.indices[3 * i + 1] = patch.indices[3 * i + 2];
            patch.indices[3 * i + 2] = temp;
        }
    }

    if (fair && params->smooth_iterations) {
        result = refine_patch(&patch);
        if (result != GEOMLIB_MESH_SUCCESS) goto finish;
        size_t added = patch.vertex_count - count;
        interior = (uint32_t*)malloc(added * sizeof(uint32_t));
        if (!interior) {
            result = GEOMLIB_MESH_ERROR_MEMORY;
            goto finish;
        }
        for (size_t i = 0; i < added; i++) interior[i] = (uint32_t)(count + i);
        result = geomlib_fair_filled_region(&patch, interior, added, params->smooth_iterations,
                                            params->smooth_lambda);
        if (result != GEOMLIB_MESH_SUCCESS) goto finish;
        if (mesh->normals) {
            geomlib_trimesh_compute_normals(&patch);
            if (!patch.normals) {
                result = GEOMLIB_MESH_ERROR_MEMORY;
                goto finish;
            }
        }
    }

    size_t base = mesh->vertex_count;
    size_t added = patch.vertex_count - count;
    if (!mesh->triangle_count || base > UINT32_MAX - added ||
        base + added > SIZE_MAX / sizeof(geomlib_vec3) ||
        mesh->triangle_count > SIZE_MAX / (3 * sizeof(uint32_t)) - patch.triangle_count) {
        result = GEOMLIB_MESH_ERROR_INVALID;
        goto finish;
    }
    result = geomlib_trimesh_reserve_vertices(mesh, base + added);
    if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    result = geomlib_trimesh_reserve_triangles(mesh, mesh->triangle_count + patch.triangle_count);
    if (result != GEOMLIB_MESH_SUCCESS) goto finish;
    for (size_t i = 0; i < added; i++) {
        mesh->vertices[base + i] = patch.vertices[count + i];
        if (mesh->normals) mesh->normals[base + i] = patch.normals[count + i];
        if (mesh->texcoords) mesh->texcoords[base + i] = patch.texcoords[count + i];
    }
    for (size_t i = 0; i < patch.triangle_count * 3; i++) {
        uint32_t local = patch.indices[i];
        mesh->indices[mesh->triangle_count * 3 + i] = local < count
            ? hole->boundary_vertices[local] : (uint32_t)(base + local - count);
    }
    mesh->vertex_count += added;
    mesh->triangle_count += patch.triangle_count;

finish:
    free(points);
    free(interior);
    geomlib_trimesh_free(&patch);
    return result;
}

int geomlib_fill_tangent_plane(geomlib_trimesh *mesh, const geomlib_hole *hole,
                               const geomlib_fill_params *params) {
    return fill_patch(mesh, hole, params, false, false);
}

int geomlib_fill_advancing_front(geomlib_trimesh *mesh, const geomlib_hole *hole,
                                 const geomlib_fill_params *params) {
    return fill_patch(mesh, hole, params, true, params && params->fair_interior);
}

int geomlib_fill_minimal_surface(geomlib_trimesh *mesh, const geomlib_hole *hole,
                                 const geomlib_fill_params *params) {
    return fill_patch(mesh, hole, params, false, params && params->fair_interior);
}

int geomlib_fill_natural(geomlib_trimesh *mesh, const geomlib_hole *hole,
                         const geomlib_fill_params *params) {
    if (!params) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_fill_params natural = *params;
    if (!natural.smooth_iterations) natural.smooth_iterations = 10;
    return fill_patch(mesh, hole, &natural, true, true);
}

/* ============================================================================
 * Hole Filling Main Functions
 * ============================================================================ */

int geomlib_fill_hole(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
) {
    if (!params) return GEOMLIB_MESH_ERROR_INVALID;
    switch (params->method) {
        case GEOMLIB_FILL_TANGENT_PLANE:
            return geomlib_fill_tangent_plane(mesh, hole, params);

        case GEOMLIB_FILL_ADVANCING_FRONT:
            return geomlib_fill_advancing_front(mesh, hole, params);

        case GEOMLIB_FILL_MINIMAL_SURFACE:
            return geomlib_fill_minimal_surface(mesh, hole, params);

        case GEOMLIB_FILL_NATURAL:
            return geomlib_fill_natural(mesh, hole, params);

        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }
}

int geomlib_fill_all_holes(
    geomlib_trimesh *mesh,
    const geomlib_fill_params *params
) {
    if (!mesh || !params) return GEOMLIB_MESH_ERROR_INVALID;
    size_t old_vertices = mesh->vertex_count, old_triangles = mesh->triangle_count;
    geomlib_holes holes = geomlib_detect_holes(mesh);

    /* Fill each hole */
    for (size_t i = 0; i < holes.hole_count; i++) {
        int result = geomlib_fill_hole(mesh, &holes.holes[i], params);
        if (result != GEOMLIB_MESH_SUCCESS) {
            mesh->vertex_count = old_vertices;
            mesh->triangle_count = old_triangles;
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

static int fair_region(geomlib_trimesh *mesh, const geomlib_hole *hole,
                       const uint32_t *vertices, size_t count, uint32_t iterations, float lambda) {
    if (!mesh || (count && !vertices) || !isfinite(lambda) || lambda < 0.0f || lambda > 1.0f)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!count || !iterations) return GEOMLIB_MESH_SUCCESS;

    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);
    geomlib_vec3 *positions = (geomlib_vec3*)malloc(count * sizeof(geomlib_vec3));
    bool *fixed = (bool*)calloc(mesh->vertex_count, sizeof(bool));
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (!positions || !fixed || edges.vertex_count != mesh->vertex_count) goto finish;
    result = GEOMLIB_MESH_ERROR_INVALID;
    for (size_t i = 0; hole && i < hole->vertex_count; i++) {
        uint32_t vertex = hole->boundary_vertices[i];
        if (vertex >= mesh->vertex_count) goto finish;
        fixed[vertex] = true;
    }
    for (size_t i = 0; i < count; i++) {
        if (vertices[i] >= mesh->vertex_count) goto finish;
    }
    for (uint32_t step = 0; step < iterations; step++) {
        for (size_t i = 0; i < count; i++) {
            uint32_t vertex = vertices[i];
            positions[i] = mesh->vertices[vertex];
            uint32_t degree = edges.vertex_edge_count[vertex];
            if (fixed[vertex] || !degree) continue;
            geomlib_vec3 sum = {0};
            for (uint32_t j = 0; j < degree; j++) {
                geomlib_edge edge = edges.edges[edges.vertex_edges[vertex][j]];
                uint32_t neighbor = edge.v0 == vertex ? edge.v1 : edge.v0;
                sum = geomlib_vec3_add(sum, mesh->vertices[neighbor]);
            }
            geomlib_vec3 delta = geomlib_vec3_sub(geomlib_vec3_mul(sum, 1.0f / degree), positions[i]);
            positions[i] = geomlib_vec3_add(positions[i], geomlib_vec3_mul(delta, lambda));
        }
        for (size_t i = 0; i < count; i++) mesh->vertices[vertices[i]] = positions[i];
    }
    result = GEOMLIB_MESH_SUCCESS;

finish:
    geomlib_edge_mesh_free(&edges);
    free(positions);
    free(fixed);
    return result;
}

int geomlib_smooth_hole_interior(geomlib_trimesh *mesh, const geomlib_hole *hole,
                                  const uint32_t *interior_vertices, size_t interior_count,
                                  uint32_t iterations) {
    if (!hole || (hole->vertex_count && !hole->boundary_vertices)) return GEOMLIB_MESH_ERROR_INVALID;
    return fair_region(mesh, hole, interior_vertices, interior_count, iterations, 0.5f);
}

int geomlib_fair_filled_region(geomlib_trimesh *mesh, const uint32_t *region_vertices,
                                size_t vertex_count, uint32_t iterations, float lambda) {
    return fair_region(mesh, NULL, region_vertices, vertex_count, iterations, lambda);
}
