/*
 * geomlib - Geodesic Distance Computation Implementation
 * Heat method with matrix-free conjugate gradient solves
 */

#include "geomlib/geomlib_geodesic.h"
#include "geomlib/geomlib_healing.h"
#include "geomlib/geomlib_curvature.h"
#include "geomlib/geomlib_solver.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* Create geodesic distance structure */
geomlib_geodesic_distances geomlib_geodesic_create(size_t vertex_count) {
    geomlib_geodesic_distances geo;
    geo.distances = calloc(vertex_count, sizeof(float));
    geo.vertex_count = vertex_count;
    geo.source_vertices = NULL;
    geo.source_count = 0;
    return geo;
}

/* Free geodesic distance structure */
void geomlib_geodesic_free(geomlib_geodesic_distances *geodesic) {
    if (!geodesic) return;
    if (geodesic->distances) {
        free(geodesic->distances);
        geodesic->distances = NULL;
    }
    if (geodesic->source_vertices) {
        free(geodesic->source_vertices);
        geodesic->source_vertices = NULL;
    }
    geodesic->vertex_count = 0;
    geodesic->source_count = 0;
}

/* Compute mean edge length */
float geomlib_compute_mean_edge_length_for_geodesic(const geomlib_trimesh *mesh) {
    if (!mesh || mesh->triangle_count == 0) return 1.0f;

    float total_length = 0.0f;
    size_t edge_count = 0;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        total_length += geomlib_vec3_length(geomlib_vec3_sub(p1, p0));
        total_length += geomlib_vec3_length(geomlib_vec3_sub(p2, p1));
        total_length += geomlib_vec3_length(geomlib_vec3_sub(p0, p2));
        edge_count += 3;
    }

    return edge_count > 0 ? (total_length / edge_count) : 1.0f;
}

/* ============================================================================
 * Heat method (Crane et al. 2012) with matrix-free conjugate gradients
 * ============================================================================ */

typedef struct geo_operator {
    size_t count;
    size_t *offsets;
    uint32_t *neighbors;
    float *weights;
    double *diag;
    double *mass;
} geo_operator;

typedef struct geo_apply {
    const geo_operator *op;
    double time_step;
} geo_apply;

typedef struct raw_geo_edge {
    uint32_t a, b;
    double weight;
} raw_geo_edge;

static int compare_raw_geo_edges(const void *left, const void *right) {
    const raw_geo_edge *x = left, *y = right;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    if (x->b != y->b) return x->b < y->b ? -1 : 1;
    return 0;
}

static void geo_operator_free(geo_operator *op) {
    if (!op) return;
    free(op->offsets);
    free(op->neighbors);
    free(op->weights);
    free(op->diag);
    free(op->mass);
    memset(op, 0, sizeof(*op));
}

/* Cotangent stiffness and lumped mass built from the triangles. */
static int geo_build_operator(const geomlib_trimesh *mesh, geo_operator *op) {
    size_t n = mesh->vertex_count, triangles = mesh->triangle_count;
    memset(op, 0, sizeof(*op));
    if (!n || !triangles || triangles > SIZE_MAX / (3 * sizeof(raw_geo_edge)))
        return GEOMLIB_MESH_ERROR_INVALID;
    raw_geo_edge *raw = malloc(triangles * 3 * sizeof(*raw));
    double *mass = calloc(n, sizeof(*mass));
    double *diag = calloc(n, sizeof(*diag));
    if (!raw || !mass || !diag) {
        free(raw);
        free(mass);
        free(diag);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    size_t raw_count = 0;
    for (size_t t = 0; t < triangles; t++) {
        uint32_t index[3] = {mesh->indices[t * 3], mesh->indices[t * 3 + 1],
                             mesh->indices[t * 3 + 2]};
        geomlib_vec3d p[3];
        for (int i = 0; i < 3; i++) {
            if (index[i] >= n) { free(raw); free(mass); free(diag);
                                 return GEOMLIB_MESH_ERROR_INVALID; }
            p[i] = (geomlib_vec3d){mesh->vertices[index[i]].x, mesh->vertices[index[i]].y,
                                   mesh->vertices[index[i]].z};
        }
        geomlib_vec3d normal = geomlib_vec3d_cross(geomlib_vec3d_sub(p[1], p[0]),
                                                   geomlib_vec3d_sub(p[2], p[0]));
        double area2 = geomlib_vec3d_length(normal);
        if (!(area2 > 1e-20)) continue;
        double area = 0.5 * area2;
        for (int i = 0; i < 3; i++) mass[index[i]] += area / 3.0;
        for (int i = 0; i < 3; i++) {
            int j = (i + 1) % 3, k = (i + 2) % 3;
            geomlib_vec3d u = geomlib_vec3d_sub(p[k], p[i]);
            geomlib_vec3d v = geomlib_vec3d_sub(p[k], p[j]);
            double cross = geomlib_vec3d_length(geomlib_vec3d_cross(u, v));
            if (!(cross > 1e-20)) continue;
            double cot = geomlib_vec3d_dot(u, v) / cross;
            raw_geo_edge entry;
            if (index[i] < index[j]) { entry.a = index[i]; entry.b = index[j]; }
            else { entry.a = index[j]; entry.b = index[i]; }
            entry.weight = 0.5 * cot;
            raw[raw_count++] = entry;
        }
    }
    qsort(raw, raw_count, sizeof(*raw), compare_raw_geo_edges);
    size_t unique = 0;
    for (size_t start = 0; start < raw_count;) {
        size_t end = start + 1;
        double weight = raw[start].weight;
        while (end < raw_count && raw[end].a == raw[start].a && raw[end].b == raw[start].b) {
            weight += raw[end].weight;
            end++;
        }
        raw[unique].a = raw[start].a;
        raw[unique].b = raw[start].b;
        raw[unique].weight = weight;
        unique++;
        start = end;
    }

    size_t *counts = calloc(n, sizeof(*counts));
    size_t *offsets = calloc(n + 1, sizeof(*offsets));
    if (!counts || !offsets) {
        free(raw); free(counts); free(offsets); free(mass); free(diag);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t e = 0; e < unique; e++) {
        counts[raw[e].a]++;
        counts[raw[e].b]++;
    }
    for (size_t i = 0; i < n; i++) offsets[i + 1] = offsets[i] + counts[i];
    size_t entries = offsets[n];
    uint32_t *neighbors = malloc(entries * sizeof(*neighbors));
    float *weights = malloc(entries * sizeof(*weights));
    size_t *cursor = malloc(n * sizeof(*cursor));
    if (!neighbors || !weights || !cursor) {
        free(raw); free(counts); free(offsets); free(neighbors); free(weights); free(cursor);
        free(mass); free(diag);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memcpy(cursor, offsets, n * sizeof(*cursor));
    for (size_t e = 0; e < unique; e++) {
        uint32_t a = raw[e].a, b = raw[e].b;
        float weight = (float)raw[e].weight;
        neighbors[cursor[a]] = b;
        weights[cursor[a]++] = weight;
        neighbors[cursor[b]] = a;
        weights[cursor[b]++] = weight;
        diag[a] += raw[e].weight;
        diag[b] += raw[e].weight;
    }
    free(raw);
    free(counts);
    free(cursor);

    op->count = n;
    op->offsets = offsets;
    op->neighbors = neighbors;
    op->weights = weights;
    op->diag = diag;
    op->mass = mass;
    return GEOMLIB_MESH_SUCCESS;
}

static void apply_heat(const void *context, const float *x, float *y) {
    const geo_apply *apply = context;
    const geo_operator *op = apply->op;
    double t = apply->time_step;
    for (size_t i = 0; i < op->count; i++) {
        double value = (op->mass[i] + t * op->diag[i]) * x[i];
        for (size_t k = op->offsets[i]; k < op->offsets[i + 1]; k++)
            value -= t * op->weights[k] * x[op->neighbors[k]];
        if (!(op->mass[i] + t * op->diag[i] > 0.0)) value = x[i];
        y[i] = (float)value;
    }
}

static void apply_poisson(const void *context, const float *x, float *y) {
    const geo_apply *apply = context;
    const geo_operator *op = apply->op;
    double epsilon = apply->time_step;
    for (size_t i = 0; i < op->count; i++) {
        double value = (op->diag[i] + epsilon) * x[i];
        for (size_t k = op->offsets[i]; k < op->offsets[i + 1]; k++)
            value -= op->weights[k] * x[op->neighbors[k]];
        y[i] = (float)value;
    }
}

typedef void (*geo_apply_fn)(const void *context, const float *x, float *y);

static double geo_dot(const float *a, const float *b, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += (double)a[i] * b[i];
    return sum;
}

static int geo_cg(size_t n, geo_apply_fn apply, const void *context, const float *b, float *x,
                  int max_iterations, float tolerance) {
    float *residual = malloc(n * sizeof(*residual));
    float *direction = malloc(n * sizeof(*direction));
    float *product = malloc(n * sizeof(*product));
    if (!residual || !direction || !product) {
        free(residual); free(direction); free(product);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    apply(context, x, product);
    for (size_t i = 0; i < n; i++) {
        residual[i] = b[i] - product[i];
        direction[i] = residual[i];
    }
    double rs = geo_dot(residual, residual, n);
    double rs0 = rs;
    int status = GEOMLIB_MESH_SUCCESS;
    if (rs0 <= 0.0) goto finish;
    for (int iter = 0; iter < max_iterations; iter++) {
        apply(context, direction, product);
        double denominator = geo_dot(direction, product, n);
        if (!(fabs(denominator) > 1e-30)) { status = GEOMLIB_SOLVER_ERROR_SINGULAR; break; }
        double alpha = rs / denominator;
        for (size_t i = 0; i < n; i++) x[i] += (float)(alpha * direction[i]);
        for (size_t i = 0; i < n; i++) residual[i] -= (float)(alpha * product[i]);
        double rs_new = geo_dot(residual, residual, n);
        if (sqrt(rs_new) <= (double)tolerance * sqrt(rs0)) { rs = rs_new; goto finish; }
        double beta = rs_new / rs;
        for (size_t i = 0; i < n; i++) direction[i] = residual[i] + (float)(beta * direction[i]);
        rs = rs_new;
    }
    status = GEOMLIB_SOLVER_ERROR_NO_CONVERGE;

finish:
    free(residual);
    free(direction);
    free(product);
    return status;
}

/* Face gradients of barycentric basis functions: grad b_i = (n x opposite_edge) / (2A). */
static void face_basis_gradient(const geomlib_trimesh *mesh, size_t triangle,
                                geomlib_vec3 grad[3], geomlib_vec3 *normal, double *area) {
    uint32_t i0 = mesh->indices[triangle * 3], i1 = mesh->indices[triangle * 3 + 1],
             i2 = mesh->indices[triangle * 3 + 2];
    geomlib_vec3 p0 = mesh->vertices[i0], p1 = mesh->vertices[i1], p2 = mesh->vertices[i2];
    geomlib_vec3 n = geomlib_vec3_cross(geomlib_vec3_sub(p1, p0), geomlib_vec3_sub(p2, p0));
    float length = geomlib_vec3_length(n);
    *area = 0.5 * (double)length;
    if (!(length > 1e-20f)) {
        grad[0] = grad[1] = grad[2] = (geomlib_vec3){0, 0, 0};
        *normal = (geomlib_vec3){0, 0, 1};
        return;
    }
    *normal = geomlib_vec3_mul(n, 1.0f / length);
    float scale = 1.0f / length;
    grad[0] = geomlib_vec3_mul(geomlib_vec3_cross(*normal, geomlib_vec3_sub(p2, p1)), scale);
    grad[1] = geomlib_vec3_mul(geomlib_vec3_cross(*normal, geomlib_vec3_sub(p0, p2)), scale);
    grad[2] = geomlib_vec3_mul(geomlib_vec3_cross(*normal, geomlib_vec3_sub(p1, p0)), scale);
}

static int geodesic_heat_solve(const geomlib_trimesh *mesh, const uint32_t *sources,
                               size_t source_count, float time_step,
                               geomlib_geodesic_distances *distances) {
    size_t n = mesh->vertex_count;
    geo_operator op;
    int status = geo_build_operator(mesh, &op);
    if (status != GEOMLIB_MESH_SUCCESS) return status;

    float *u0 = calloc(n, sizeof(*u0));
    float *u = calloc(n, sizeof(*u));
    float *b = calloc(n, sizeof(*b));
    float *divergence = calloc(n, sizeof(*divergence));
    float *phi = calloc(n, sizeof(*phi));
    geomlib_vec3 *face_vectors = malloc(mesh->triangle_count * sizeof(*face_vectors));
    if (!u0 || !u || !b || !divergence || !phi || !face_vectors) {
        status = GEOMLIB_MESH_ERROR_MEMORY;
        goto finish;
    }
    for (size_t i = 0; i < source_count; i++)
        if (sources[i] < n) u0[sources[i]] = 1.0f;
    for (size_t i = 0; i < n; i++) {
        u[i] = u0[i];
        b[i] = (float)(op.mass[i] * (double)u0[i]);
    }

    double mean_edge = geomlib_compute_mean_edge_length_for_geodesic(mesh);
    double t = time_step > 0.0f && isfinite(time_step) ? time_step : mean_edge * mean_edge;
    if (!(t > 0.0)) t = 1e-6;
    geo_apply heat = {&op, t};
    status = geo_cg(n, apply_heat, &heat, b, u, 500, 1e-6f);
    if (status != GEOMLIB_MESH_SUCCESS) goto finish;

    for (size_t tri = 0; tri < mesh->triangle_count; tri++) {
        geomlib_vec3 grad[3], normal;
        double area;
        face_basis_gradient(mesh, tri, grad, &normal, &area);
        uint32_t i0 = mesh->indices[tri * 3], i1 = mesh->indices[tri * 3 + 1],
                 i2 = mesh->indices[tri * 3 + 2];
        geomlib_vec3 grad_u = geomlib_vec3_add(
            geomlib_vec3_add(geomlib_vec3_mul(grad[0], u[i0]),
                             geomlib_vec3_mul(grad[1], u[i1])),
            geomlib_vec3_mul(grad[2], u[i2]));
        float length = geomlib_vec3_length(grad_u);
        face_vectors[tri] = length > 1e-20f
            ? geomlib_vec3_mul(grad_u, -1.0f / length) : (geomlib_vec3){0, 0, 0};
        if (!(area > 0.0)) continue;
        divergence[i0] += (float)(area * geomlib_vec3_dot(grad[0], face_vectors[tri]));
        divergence[i1] += (float)(area * geomlib_vec3_dot(grad[1], face_vectors[tri]));
        divergence[i2] += (float)(area * geomlib_vec3_dot(grad[2], face_vectors[tri]));
    }

    double trace = 0.0;
    for (size_t i = 0; i < n; i++) trace += op.diag[i];
    double epsilon = 1e-8 * (n ? trace / (double)n : 1.0) + 1e-12;
    geo_apply poisson = {&op, epsilon};
    status = geo_cg(n, apply_poisson, &poisson, divergence, phi, 500, 1e-6f);
    if (status != GEOMLIB_MESH_SUCCESS) goto finish;

    double source_phi = 0.0;
    size_t valid_sources = 0;
    for (size_t i = 0; i < source_count; i++)
        if (sources[i] < n) { source_phi += phi[sources[i]]; valid_sources++; }
    if (valid_sources) source_phi /= (double)valid_sources;
    for (size_t i = 0; i < n; i++) {
        float value = (float)((double)phi[i] - source_phi);
        if (value < 0.0f) value = 0.0f;
        distances->distances[i] = value;
    }
    status = GEOMLIB_MESH_SUCCESS;

finish:
    geo_operator_free(&op);
    free(u0);
    free(u);
    free(b);
    free(divergence);
    free(phi);
    free(face_vectors);
    return status;
}

/* Edge-following fallback, used only if the linear solves fail. */
static int geodesic_edge_propagation(const geomlib_trimesh *mesh, const uint32_t *sources,
                                     size_t source_count, geomlib_geodesic_distances *distances) {
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;
    for (size_t i = 0; i < mesh->vertex_count; i++) distances->distances[i] = FLT_MAX;
    for (size_t i = 0; i < source_count; i++)
        if (sources[i] < mesh->vertex_count) distances->distances[sources[i]] = 0.0f;
    bool *visited = calloc(mesh->vertex_count, sizeof(bool));
    if (!visited) { geomlib_halfedge_free(&hemesh); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (size_t iter = 0; iter < mesh->vertex_count; iter++) {
        uint32_t min_vertex = UINT32_MAX;
        float min_distance = FLT_MAX;
        for (size_t v = 0; v < mesh->vertex_count; v++)
            if (!visited[v] && distances->distances[v] < min_distance) {
                min_distance = distances->distances[v];
                min_vertex = (uint32_t)v;
            }
        if (min_vertex == UINT32_MAX) break;
        visited[min_vertex] = true;
        for (size_t he = 0; he < hemesh.halfedge_count; he++) {
            if (hemesh.halfedges[he].vertex != min_vertex) continue;
            uint32_t next = hemesh.halfedges[he].next;
            if (next == UINT32_MAX) continue;
            uint32_t neighbor = hemesh.halfedges[next].vertex;
            if (visited[neighbor]) continue;
            float length = geomlib_vec3_length(
                geomlib_vec3_sub(mesh->vertices[neighbor], mesh->vertices[min_vertex]));
            float candidate = distances->distances[min_vertex] + length;
            if (candidate < distances->distances[neighbor])
                distances->distances[neighbor] = candidate;
        }
    }
    free(visited);
    geomlib_halfedge_free(&hemesh);
    return GEOMLIB_MESH_SUCCESS;
}

/* Geodesic distances by the heat method. time_step <= 0 selects mean_edge_length^2. */
int geomlib_compute_geodesic_heat(
    const geomlib_trimesh *mesh,
    const uint32_t *source_vertices,
    size_t source_count,
    float time_step,
    geomlib_geodesic_distances *distances
) {
    if (!mesh || !source_vertices || source_count == 0 || !distances ||
        distances->vertex_count != mesh->vertex_count)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!isfinite(time_step)) return GEOMLIB_MESH_ERROR_INVALID;

    if (distances->source_vertices) free(distances->source_vertices);
    distances->source_vertices = malloc(source_count * sizeof(uint32_t));
    if (!distances->source_vertices) return GEOMLIB_MESH_ERROR_MEMORY;
    memcpy(distances->source_vertices, source_vertices, source_count * sizeof(uint32_t));
    distances->source_count = source_count;

    int status = geodesic_heat_solve(mesh, source_vertices, source_count, time_step, distances);
    if (status != GEOMLIB_MESH_SUCCESS)
        status = geodesic_edge_propagation(mesh, source_vertices, source_count, distances);
    return status;
}

/* Compute geodesic distance from a single source vertex */
int geomlib_compute_geodesic_from_vertex(
    const geomlib_trimesh *mesh,
    uint32_t source_vertex,
    geomlib_geodesic_distances *distances
) {
    if (source_vertex >= mesh->vertex_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    return geomlib_compute_geodesic_heat(mesh, &source_vertex, 1, 0.0f, distances);
}

/* Compute all-pairs geodesic distances */
int geomlib_compute_geodesic_all_pairs(
    const geomlib_trimesh *mesh,
    float **distance_matrix
) {
    if (!mesh || !distance_matrix) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Allocate distance matrix */
    size_t n = mesh->vertex_count;
    *distance_matrix = malloc(n * n * sizeof(float));
    if (!*distance_matrix) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Compute distances from each vertex */
    geomlib_geodesic_distances geo = geomlib_geodesic_create(n);

    for (size_t i = 0; i < n; i++) {
        uint32_t source = (uint32_t)i;
        int result = geomlib_compute_geodesic_from_vertex(mesh, source, &geo);

        if (result != GEOMLIB_MESH_SUCCESS) {
            free(*distance_matrix);
            *distance_matrix = NULL;
            geomlib_geodesic_free(&geo);
            return result;
        }

        /* Copy distances to matrix row */
        for (size_t j = 0; j < n; j++) {
            (*distance_matrix)[i * n + j] = geo.distances[j];
        }
    }

    geomlib_geodesic_free(&geo);

    return GEOMLIB_MESH_SUCCESS;
}
