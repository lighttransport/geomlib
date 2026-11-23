/*
 * geomlib - QEM Decimation Implementation
 * Quadric-error edge-collapse decimation (Garland-Heckbert)
 */

#include "geomlib/geomlib_decimation.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Create default decimation parameters */
geomlib_decimation_params geomlib_decimation_params_default(float target_reduction) {
    geomlib_decimation_params params;
    params.target_reduction = target_reduction;
    params.target_vertex_count = 0;
    params.target_triangle_count = 0;
    params.max_error = 0.0f;  /* No limit */
    params.preserve_boundaries = true;
    params.preserve_topology = true;
    params.boundary_weight = 1000.0f;
    return params;
}

/* Initialize a quadric from a plane */
geomlib_quadric geomlib_quadric_from_plane(float a, float b, float c, float d) {
    geomlib_quadric q;
    q.a2 = a * a;
    q.ab = a * b;
    q.ac = a * c;
    q.ad = a * d;
    q.b2 = b * b;
    q.bc = b * c;
    q.bd = b * d;
    q.c2 = c * c;
    q.cd = c * d;
    q.d2 = d * d;
    return q;
}

/* Add two quadrics */
geomlib_quadric geomlib_quadric_add(geomlib_quadric q1, geomlib_quadric q2) {
    geomlib_quadric result;
    result.a2 = q1.a2 + q2.a2;
    result.ab = q1.ab + q2.ab;
    result.ac = q1.ac + q2.ac;
    result.ad = q1.ad + q2.ad;
    result.b2 = q1.b2 + q2.b2;
    result.bc = q1.bc + q2.bc;
    result.bd = q1.bd + q2.bd;
    result.c2 = q1.c2 + q2.c2;
    result.cd = q1.cd + q2.cd;
    result.d2 = q1.d2 + q2.d2;
    return result;
}

/* Compute quadric error for a vertex position */
double geomlib_quadric_error(geomlib_quadric q, geomlib_vec3 v) {
    double x = v.x, y = v.y, z = v.z;

    double error =
        q.a2 * x * x + 2.0 * q.ab * x * y + 2.0 * q.ac * x * z + 2.0 * q.ad * x +
        q.b2 * y * y + 2.0 * q.bc * y * z + 2.0 * q.bd * y +
        q.c2 * z * z + 2.0 * q.cd * z +
        q.d2;

    return error;
}

/* ============================================================================
 * Quadric error edge-collapse decimation
 * ============================================================================ */

typedef struct qem_raw_edge {
    uint32_t a, b;
    uint32_t opposite;
} qem_raw_edge;

typedef struct qem_edge {
    uint32_t a, b;
    uint32_t opposite[2];
    uint8_t faces;
    double cost;
    geomlib_vec3 target;
} qem_edge;

static int compare_qem_raw_edges(const void *left, const void *right) {
    const qem_raw_edge *x = left, *y = right;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    if (x->b != y->b) return x->b < y->b ? -1 : 1;
    return (x->opposite > y->opposite) - (x->opposite < y->opposite);
}

static const qem_edge *qem_find_edge(const qem_edge *edges, size_t count,
                                     uint32_t x, uint32_t y) {
    if (x > y) { uint32_t t = x; x = y; y = t; }
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (edges[mid].a < x || (edges[mid].a == x && edges[mid].b < y)) lo = mid + 1;
        else hi = mid;
    }
    return lo < count && edges[lo].a == x && edges[lo].b == y ? &edges[lo] : NULL;
}

/* Solve the quadric for its optimal position, falling back to the endpoints and midpoint. */
static geomlib_vec3 qem_optimal_position(geomlib_quadric q, geomlib_vec3 a, geomlib_vec3 b,
                                         double *cost) {
    double a2 = q.a2, ab = q.ab, ac = q.ac, b2 = q.b2, bc = q.bc, c2 = q.c2;
    double det = a2 * (b2 * c2 - bc * bc) - ab * (ab * c2 - bc * ac) +
                 ac * (ab * bc - b2 * ac);
    if (fabs(det) > 1e-18) {
        double r0 = -q.ad, r1 = -q.bd, r2 = -q.cd;
        double x = (r0 * (b2 * c2 - bc * bc) - ab * (r1 * c2 - bc * r2) +
                    ac * (r1 * bc - b2 * r2)) / det;
        double y = (a2 * (r1 * c2 - bc * r2) - r0 * (ab * c2 - bc * ac) +
                    ac * (ab * r2 - r1 * ac)) / det;
        double z = (a2 * (b2 * r2 - r1 * bc) - ab * (ab * r2 - r1 * ac) +
                    r0 * (ab * bc - b2 * ac)) / det;
        if (isfinite(x) && isfinite(y) && isfinite(z)) {
            geomlib_vec3 point = {(float)x, (float)y, (float)z};
            *cost = geomlib_quadric_error(q, point);
            if (isfinite(*cost)) return point;
        }
    }
    geomlib_vec3 midpoint = geomlib_vec3_mul(geomlib_vec3_add(a, b), 0.5f);
    geomlib_vec3 candidates[3] = {a, b, midpoint};
    geomlib_vec3 best = midpoint;
    double best_cost = DBL_MAX;
    for (int i = 0; i < 3; i++) {
        double value = geomlib_quadric_error(q, candidates[i]);
        if (value < best_cost) { best_cost = value; best = candidates[i]; }
    }
    *cost = best_cost;
    return best;
}

static int qem_build_edges(const geomlib_trimesh *mesh, qem_edge **edges_out, size_t *count_out) {
    size_t triangle_count = mesh->triangle_count;
    if (!triangle_count) { *edges_out = NULL; *count_out = 0; return GEOMLIB_MESH_SUCCESS; }
    if (triangle_count > SIZE_MAX / (3 * sizeof(qem_raw_edge)))
        return GEOMLIB_MESH_ERROR_MEMORY;
    qem_raw_edge *raw = malloc(triangle_count * 3 * sizeof(*raw));
    qem_edge *edges = malloc(triangle_count * 3 * sizeof(*edges));
    if (!raw || !edges) { free(raw); free(edges); return GEOMLIB_MESH_ERROR_MEMORY; }
    size_t raw_count = 0;
    for (size_t t = 0; t < triangle_count; t++) {
        uint32_t v[3] = {mesh->indices[t * 3], mesh->indices[t * 3 + 1],
                         mesh->indices[t * 3 + 2]};
        for (int j = 0; j < 3; j++) {
            uint32_t x = v[j], y = v[(j + 1) % 3], opposite = v[(j + 2) % 3];
            if (x == y) continue;
            qem_raw_edge entry;
            if (x < y) { entry.a = x; entry.b = y; } else { entry.a = y; entry.b = x; }
            entry.opposite = opposite;
            raw[raw_count++] = entry;
        }
    }
    qsort(raw, raw_count, sizeof(*raw), compare_qem_raw_edges);
    size_t unique = 0;
    for (size_t start = 0; start < raw_count;) {
        size_t end = start + 1;
        while (end < raw_count && raw[end].a == raw[start].a && raw[end].b == raw[start].b) end++;
        qem_edge edge;
        edge.a = raw[start].a;
        edge.b = raw[start].b;
        edge.opposite[0] = raw[start].opposite;
        edge.opposite[1] = UINT32_MAX;
        if (end - start > 1) edge.opposite[1] = raw[start + 1].opposite;
        edge.faces = (uint8_t)(end - start > 2 ? 3 : end - start);
        edge.cost = DBL_MAX;
        edge.target = (geomlib_vec3){0, 0, 0};
        edges[unique++] = edge;
        start = end;
    }
    free(raw);
    *edges_out = edges;
    *count_out = unique;
    return GEOMLIB_MESH_SUCCESS;
}

typedef struct qem_face {
    uint32_t key[3];
    uint32_t v[3];
} qem_face;

static int compare_qem_faces(const void *left, const void *right) {
    const qem_face *x = left, *y = right;
    for (int i = 0; i < 3; i++)
        if (x->key[i] != y->key[i]) return x->key[i] < y->key[i] ? -1 : 1;
    return 0;
}

int geomlib_decimate_qem(
    geomlib_trimesh *mesh,
    const geomlib_decimation_params *params,
    geomlib_decimation_stats *stats
) {
    if (!mesh || !params) return GEOMLIB_MESH_ERROR_INVALID;

    size_t original_vertices = mesh->vertex_count;
    size_t original_triangles = mesh->triangle_count;
    size_t target_triangles;
    if (params->target_triangle_count > 0) {
        target_triangles = params->target_triangle_count;
    } else if (params->target_vertex_count > 0) {
        target_triangles = params->target_vertex_count > 2
            ? (params->target_vertex_count - 2) * 2 : 4;
    } else {
        double reduction = params->target_reduction;
        if (reduction < 0.0) reduction = 0.0;
        if (reduction > 1.0) reduction = 1.0;
        target_triangles = (size_t)((1.0 - reduction) * (double)original_triangles);
    }
    if (target_triangles < 4) target_triangles = 4;
    if (target_triangles >= original_triangles) {
        if (stats) {
            memset(stats, 0, sizeof(*stats));
            stats->original_vertices = original_vertices;
            stats->original_triangles = original_triangles;
            stats->final_vertices = original_vertices;
            stats->final_triangles = original_triangles;
        }
        return GEOMLIB_MESH_SUCCESS;
    }

    geomlib_decimation_stats local_stats;
    memset(&local_stats, 0, sizeof(local_stats));
    local_stats.original_vertices = original_vertices;
    local_stats.original_triangles = original_triangles;

    geomlib_quadric *quadrics = calloc(mesh->vertex_count, sizeof(*quadrics));
    bool *removed = calloc(mesh->vertex_count, sizeof(*removed));
    if (!quadrics || !removed) {
        free(quadrics); free(removed);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3], i1 = mesh->indices[t * 3 + 1],
                 i2 = mesh->indices[t * 3 + 2];
        geomlib_vec3 p0 = mesh->vertices[i0], p1 = mesh->vertices[i1], p2 = mesh->vertices[i2];
        geomlib_vec3 normal = geomlib_vec3_normalize(geomlib_vec3_cross(
            geomlib_vec3_sub(p1, p0), geomlib_vec3_sub(p2, p0)));
        float d = -(normal.x * p0.x + normal.y * p0.y + normal.z * p0.z);
        geomlib_quadric plane = geomlib_quadric_from_plane(normal.x, normal.y, normal.z, d);
        quadrics[i0] = geomlib_quadric_add(quadrics[i0], plane);
        quadrics[i1] = geomlib_quadric_add(quadrics[i1], plane);
        quadrics[i2] = geomlib_quadric_add(quadrics[i2], plane);
    }

    size_t collapsed = 0;
    double error_sum = 0.0;
    size_t max_iterations = mesh->vertex_count;
    while (mesh->triangle_count > target_triangles && collapsed < max_iterations) {
        qem_edge *edges = NULL;
        size_t edge_count = 0;
        int status = qem_build_edges(mesh, &edges, &edge_count);
        if (status != GEOMLIB_MESH_SUCCESS) {
            free(quadrics); free(removed);
            return status;
        }
        bool *boundary = NULL;
        if (params->preserve_boundaries) {
            boundary = calloc(mesh->vertex_count, sizeof(bool));
            if (!boundary) { free(edges); free(quadrics); free(removed);
                             return GEOMLIB_MESH_ERROR_MEMORY; }
            for (size_t e = 0; e < edge_count; e++)
                if (edges[e].faces == 1) {
                    boundary[edges[e].a] = true;
                    boundary[edges[e].b] = true;
                }
        }

        const qem_edge *best = NULL;
        geomlib_vec3 best_target = {0, 0, 0};
        double best_cost = DBL_MAX;
        for (size_t e = 0; e < edge_count; e++) {
            qem_edge *edge = &edges[e];
            if (edge->faces == 0 || edge->faces > 2) continue;
            if (removed[edge->a] || removed[edge->b]) continue;
            if (boundary && (boundary[edge->a] || boundary[edge->b])) continue;
            if (params->preserve_topology && edge->faces == 2) {
                uint32_t x = edge->opposite[0], y = edge->opposite[1];
                if (x == y || qem_find_edge(edges, edge_count, x, y)) continue;
            }
            geomlib_quadric sum = geomlib_quadric_add(quadrics[edge->a], quadrics[edge->b]);
            double cost;
            geomlib_vec3 target = qem_optimal_position(sum, mesh->vertices[edge->a],
                                                       mesh->vertices[edge->b], &cost);
            if (!isfinite(cost)) continue;
            if (cost < 0.0) cost = 0.0;
            if (params->max_error > 0.0f && cost > params->max_error) continue;
            if (cost < best_cost) {
                best_cost = cost;
                best = edge;
                best_target = target;
            }
        }
        if (!best) {
            free(boundary);
            free(edges);
            break;
        }

        uint32_t keep = best->a, drop = best->b;
        mesh->vertices[keep] = best_target;
        quadrics[keep] = geomlib_quadric_add(quadrics[keep], quadrics[drop]);
        removed[drop] = true;
        local_stats.edges_collapsed++;
        local_stats.max_error = (float)fmax(local_stats.max_error, best_cost);
        error_sum += best_cost;
        collapsed++;
        free(boundary);
        free(edges);

        /* Remap and compact triangles, dropping degenerates and duplicates. */
        size_t write = 0;
        for (size_t t = 0; t < mesh->triangle_count; t++) {
            uint32_t v[3] = {mesh->indices[t * 3], mesh->indices[t * 3 + 1],
                             mesh->indices[t * 3 + 2]};
            for (int j = 0; j < 3; j++) if (v[j] == drop) v[j] = keep;
            if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) continue;
            mesh->indices[write * 3] = v[0];
            mesh->indices[write * 3 + 1] = v[1];
            mesh->indices[write * 3 + 2] = v[2];
            write++;
        }
        mesh->triangle_count = write;

        /* Drop duplicate faces created by the collapse. */
        size_t faces = mesh->triangle_count;
        qem_face *list = malloc(faces * sizeof(*list));
        if (list && faces > 1) {
            for (size_t t = 0; t < faces; t++) {
                list[t].v[0] = mesh->indices[t * 3];
                list[t].v[1] = mesh->indices[t * 3 + 1];
                list[t].v[2] = mesh->indices[t * 3 + 2];
                list[t].key[0] = list[t].v[0];
                list[t].key[1] = list[t].v[1];
                list[t].key[2] = list[t].v[2];
                for (int i = 0; i < 3; i++)
                    for (int j = i + 1; j < 3; j++)
                        if (list[t].key[j] < list[t].key[i]) {
                            uint32_t swap = list[t].key[i];
                            list[t].key[i] = list[t].key[j];
                            list[t].key[j] = swap;
                        }
            }
            qsort(list, faces, sizeof(*list), compare_qem_faces);
            size_t unique = 0;
            for (size_t t = 0; t < faces; t++) {
                if (unique && list[unique - 1].key[0] == list[t].key[0] &&
                    list[unique - 1].key[1] == list[t].key[1] &&
                    list[unique - 1].key[2] == list[t].key[2]) continue;
                list[unique++] = list[t];
            }
            for (size_t t = 0; t < unique; t++) {
                mesh->indices[t * 3] = list[t].v[0];
                mesh->indices[t * 3 + 1] = list[t].v[1];
                mesh->indices[t * 3 + 2] = list[t].v[2];
            }
            mesh->triangle_count = unique;
        }
        free(list);
    }

    if (collapsed) local_stats.mean_error = (float)(error_sum / (double)collapsed);

    /* Compact removed vertices. */
    uint32_t *vertex_map = malloc(mesh->vertex_count * sizeof(*vertex_map));
    if (!vertex_map) { free(quadrics); free(removed); return GEOMLIB_MESH_ERROR_MEMORY; }
    uint32_t new_vertex_count = 0;
    for (size_t v = 0; v < mesh->vertex_count; v++)
        vertex_map[v] = removed[v] ? UINT32_MAX : new_vertex_count++;
    geomlib_vec3 *new_vertices = malloc((new_vertex_count ? new_vertex_count : 1) *
                                        sizeof(*new_vertices));
    if (!new_vertices) { free(vertex_map); free(quadrics); free(removed);
                         return GEOMLIB_MESH_ERROR_MEMORY; }
    for (size_t v = 0; v < mesh->vertex_count; v++)
        if (!removed[v]) new_vertices[vertex_map[v]] = mesh->vertices[v];
    size_t new_triangle_count = 0;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3], i1 = mesh->indices[t * 3 + 1],
                 i2 = mesh->indices[t * 3 + 2];
        if (i0 >= mesh->vertex_count || i1 >= mesh->vertex_count || i2 >= mesh->vertex_count)
            continue;
        if (vertex_map[i0] == UINT32_MAX || vertex_map[i1] == UINT32_MAX ||
            vertex_map[i2] == UINT32_MAX) continue;
        mesh->indices[new_triangle_count * 3] = vertex_map[i0];
        mesh->indices[new_triangle_count * 3 + 1] = vertex_map[i1];
        mesh->indices[new_triangle_count * 3 + 2] = vertex_map[i2];
        new_triangle_count++;
    }

    free(mesh->vertices);
    if (mesh->normals) { free(mesh->normals); mesh->normals = NULL; }
    if (mesh->texcoords) { free(mesh->texcoords); mesh->texcoords = NULL; }
    mesh->vertices = new_vertices;
    mesh->vertex_count = new_vertex_count;
    mesh->vertex_capacity = new_vertex_count;
    mesh->triangle_count = new_triangle_count;
    mesh->triangle_capacity = new_triangle_count;

    local_stats.final_vertices = new_vertex_count;
    local_stats.final_triangles = new_triangle_count;
    local_stats.reduction_ratio = original_vertices
        ? 1.0f - (float)new_vertex_count / (float)original_vertices : 0.0f;
    if (stats) *stats = local_stats;

    free(vertex_map);
    free(quadrics);
    free(removed);
    return GEOMLIB_MESH_SUCCESS;
}

/* Simplified interfaces */
int geomlib_decimate_to_ratio(geomlib_trimesh *mesh, float target_ratio) {
    geomlib_decimation_params params = geomlib_decimation_params_default(1.0f - target_ratio);
    return geomlib_decimate_qem(mesh, &params, NULL);
}

int geomlib_decimate_to_vertex_count(geomlib_trimesh *mesh, size_t target_vertices) {
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_vertex_count = target_vertices;
    return geomlib_decimate_qem(mesh, &params, NULL);
}

int geomlib_decimate_to_triangle_count(geomlib_trimesh *mesh, size_t target_triangles) {
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_triangle_count = target_triangles;
    return geomlib_decimate_qem(mesh, &params, NULL);
}
