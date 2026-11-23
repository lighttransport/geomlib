/*
 * geomlib - CFD Utilities Implementation
 * Gradient computation, interpolation, streamlines for CFD
 */

#include "geomlib/geomlib_cfd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Field Management
 * ============================================================================ */

geomlib_scalar_field geomlib_scalar_field_create(size_t count) {
    geomlib_scalar_field field;
    field.values = calloc(count, sizeof(float));
    field.count = count;
    field.gradients = NULL;
    return field;
}

void geomlib_scalar_field_free(geomlib_scalar_field *field) {
    if (!field) return;
    free(field->values);
    free(field->gradients);
    memset(field, 0, sizeof(geomlib_scalar_field));
}

geomlib_vector_field geomlib_vector_field_create(size_t count) {
    geomlib_vector_field field;
    field.values = calloc(count, sizeof(geomlib_vec3));
    field.count = count;
    return field;
}

void geomlib_vector_field_free(geomlib_vector_field *field) {
    if (!field) return;
    free(field->values);
    memset(field, 0, sizeof(geomlib_vector_field));
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

geomlib_vec3 geomlib_tetrahedron_centroid(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
) {
    geomlib_vec3 centroid;
    centroid.x = (v0.x + v1.x + v2.x + v3.x) * 0.25f;
    centroid.y = (v0.y + v1.y + v2.y + v3.y) * 0.25f;
    centroid.z = (v0.z + v1.z + v2.z + v3.z) * 0.25f;
    return centroid;
}

void geomlib_face_area_normal(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float *area,
    geomlib_vec3 *normal
) {
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);

    float len = geomlib_vec3_length(cross);
    *area = len * 0.5f;

    if (len > 1e-10f) {
        *normal = geomlib_vec3_normalize(cross);
    } else {
        *normal = geomlib_vec3_create(0.0f, 1.0f, 0.0f);
    }
}

float geomlib_point_plane_distance(
    geomlib_vec3 point,
    geomlib_vec3 plane_point,
    geomlib_vec3 plane_normal
) {
    geomlib_vec3 diff = geomlib_vec3_sub(point, plane_point);
    return geomlib_vec3_dot(diff, plane_normal);
}

/* ============================================================================
 * Gradient Computation
 * ============================================================================ */

/* Green-Gauss gradient for a polygon region */
geomlib_vec3 geomlib_gradient_green_gauss(
    const geomlib_vec3 *vertices,
    const float *values,
    size_t count
) {
    if (!vertices || !values || count < 3) return (geomlib_vec3){0};
    geomlib_vec3d normal = {0}, circulation = {0};
    geomlib_vec3d origin = {vertices[0].x, vertices[0].y, vertices[0].z};
    for (size_t i = 0; i < count; i++) {
        size_t j = (i + 1) % count;
        geomlib_vec3d a = {vertices[i].x, vertices[i].y, vertices[i].z};
        geomlib_vec3d b = {vertices[j].x, vertices[j].y, vertices[j].z};
        normal = geomlib_vec3d_add(normal, geomlib_vec3d_cross(
            geomlib_vec3d_sub(a, origin), geomlib_vec3d_sub(b, origin)));
        double mean = ((double)values[i] + values[j]) * 0.5 - values[0];
        circulation = geomlib_vec3d_add(circulation,
            geomlib_vec3d_mul(geomlib_vec3d_sub(b, a), mean));
    }
    double area2 = geomlib_vec3d_length(normal);
    if (!(area2 > 0)) return (geomlib_vec3){0};
    geomlib_vec3d gradient = geomlib_vec3d_mul(geomlib_vec3d_cross(circulation, normal),
                                              2.0 / (area2 * area2));
    return (geomlib_vec3){(float)gradient.x, (float)gradient.y, (float)gradient.z};
}

/* Least-squares gradient reconstruction */
geomlib_vec3 geomlib_gradient_least_squares(
    geomlib_vec3 center,
    const geomlib_vec3 *neighbor_positions,
    const float *neighbor_values,
    size_t neighbor_count
) {
    if (neighbor_count < 3) {
        return geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    }

    /* Build least-squares matrix: A^T * A * grad = A^T * b */
    /* where A is [dx dy dz] for each neighbor */
    float ata[3][3] = {{0}};  /* A^T * A */
    float atb[3] = {0};       /* A^T * b */

    for (size_t i = 0; i < neighbor_count; i++) {
        geomlib_vec3 delta = geomlib_vec3_sub(neighbor_positions[i], center);
        float value_diff = neighbor_values[i];

        /* Weight by inverse distance squared */
        float dist = geomlib_vec3_length(delta);
        float weight = (dist > 1e-10f) ? (1.0f / (dist * dist)) : 0.0f;

        /* Build A^T * A */
        ata[0][0] += weight * delta.x * delta.x;
        ata[0][1] += weight * delta.x * delta.y;
        ata[0][2] += weight * delta.x * delta.z;
        ata[1][1] += weight * delta.y * delta.y;
        ata[1][2] += weight * delta.y * delta.z;
        ata[2][2] += weight * delta.z * delta.z;

        /* Build A^T * b */
        atb[0] += weight * delta.x * value_diff;
        atb[1] += weight * delta.y * value_diff;
        atb[2] += weight * delta.z * value_diff;
    }

    /* Symmetric matrix */
    ata[1][0] = ata[0][1];
    ata[2][0] = ata[0][2];
    ata[2][1] = ata[1][2];

    /* Solve 3x3 system using Cramer's rule (simplified) */
    float det = ata[0][0] * (ata[1][1] * ata[2][2] - ata[1][2] * ata[2][1])
              - ata[0][1] * (ata[1][0] * ata[2][2] - ata[1][2] * ata[2][0])
              + ata[0][2] * (ata[1][0] * ata[2][1] - ata[1][1] * ata[2][0]);

    if (fabsf(det) < 1e-10f) {
        return geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    }

    geomlib_vec3 gradient;

    /* Cramer's rule for x */
    gradient.x = (atb[0] * (ata[1][1] * ata[2][2] - ata[1][2] * ata[2][1])
                - atb[1] * (ata[0][1] * ata[2][2] - ata[0][2] * ata[2][1])
                + atb[2] * (ata[0][1] * ata[1][2] - ata[0][2] * ata[1][1])) / det;

    /* Cramer's rule for y */
    gradient.y = (ata[0][0] * (atb[1] * ata[2][2] - atb[2] * ata[1][2])
                - ata[0][1] * (atb[0] * ata[2][2] - atb[2] * ata[0][2])
                + ata[0][2] * (atb[0] * ata[1][2] - atb[1] * ata[0][2])) / det;

    /* Cramer's rule for z */
    gradient.z = (ata[0][0] * (ata[1][1] * atb[2] - ata[1][2] * atb[1])
                - ata[0][1] * (ata[1][0] * atb[2] - ata[1][2] * atb[0])
                + ata[0][2] * (ata[1][0] * atb[1] - ata[1][1] * atb[0])) / det;

    return gradient;
}

static geomlib_vec3d to_double(geomlib_vec3 v) {
    return (geomlib_vec3d){v.x, v.y, v.z};
}

static geomlib_vec3 to_float(geomlib_vec3d v) {
    return (geomlib_vec3){(float)v.x, (float)v.y, (float)v.z};
}

static bool finite_vector(geomlib_vec3 v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

static bool valid_surface(const geomlib_trimesh *mesh, size_t count) {
    if (!mesh || mesh->vertex_count != count || count > UINT32_MAX ||
        mesh->triangle_count > SIZE_MAX / 3 ||
        (count && !mesh->vertices) || (mesh->triangle_count && !mesh->indices)) return false;
    for (size_t i = 0; i < count; i++) {
        if (!finite_vector(mesh->vertices[i])) return false;
    }
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        if (mesh->indices[i] >= count) return false;
    }
    return true;
}

typedef struct surface_triangle {
    geomlib_vec3d gradients[3];
    geomlib_vec3d normal;
    double area;
} surface_triangle;

static bool triangle_basis(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c,
                           surface_triangle *triangle) {
    geomlib_vec3d edge[2] = {geomlib_vec3d_sub(to_double(b), to_double(a)),
                             geomlib_vec3d_sub(to_double(c), to_double(a))};
    geomlib_vec3d cross = geomlib_vec3d_cross(edge[0], edge[1]);
    double length = geomlib_vec3d_length(cross);
    if (!(length > DBL_EPSILON * geomlib_vec3d_length(edge[0]) * geomlib_vec3d_length(edge[1])))
        return false;
    triangle->normal = geomlib_vec3d_mul(cross, 1.0 / length);
    triangle->area = 0.5 * length;
    triangle->gradients[1] = geomlib_vec3d_mul(
        geomlib_vec3d_cross(edge[1], triangle->normal), 1.0 / length);
    triangle->gradients[2] = geomlib_vec3d_mul(
        geomlib_vec3d_cross(triangle->normal, edge[0]), 1.0 / length);
    triangle->gradients[0] = geomlib_vec3d_mul(
        geomlib_vec3d_add(triangle->gradients[1], triangle->gradients[2]), -1);
    return true;
}

static bool valid_gradient_method(geomlib_gradient_method method) {
    return method == GEOMLIB_GRADIENT_GREEN_GAUSS ||
           method == GEOMLIB_GRADIENT_GREEN_GAUSS_NODE ||
           method == GEOMLIB_GRADIENT_LEAST_SQUARES;
}

static bool valid_scalar_values(const geomlib_scalar_field *field) {
    if (!field || (field->count && !field->values)) return false;
    for (size_t i = 0; i < field->count; i++) {
        if (!isfinite(field->values[i])) return false;
    }
    return true;
}

static bool valid_vector_values(const geomlib_vector_field *field) {
    if (!field || (field->count && !field->values)) return false;
    for (size_t i = 0; i < field->count; i++) {
        if (!finite_vector(field->values[i])) return false;
    }
    return true;
}

typedef enum surface_operator {
    SURFACE_GRADIENT,
    SURFACE_DIVERGENCE,
    SURFACE_CURL,
    SURFACE_LAPLACIAN
} surface_operator;

typedef struct surface_accumulator {
    double area;
    geomlib_vec3d value;
} surface_accumulator;

static int compute_surface_operator(const geomlib_trimesh *mesh, size_t count,
                                    const float *scalars, const geomlib_vec3 *vectors,
                                    surface_operator operation, float *scalar_out,
                                    geomlib_vec3 *vector_out) {
    if (!count) return GEOMLIB_MESH_SUCCESS;
    surface_accumulator *sums = calloc(count, sizeof(*sums));
    if (!sums) return GEOMLIB_MESH_ERROR_MEMORY;
    for (size_t face = 0; face < mesh->triangle_count; face++) {
        const uint32_t *indices = &mesh->indices[face * 3];
        surface_triangle triangle;
        if (!triangle_basis(mesh->vertices[indices[0]], mesh->vertices[indices[1]],
                            mesh->vertices[indices[2]], &triangle)) continue;

        geomlib_vec3d result = {0};
        if (scalars) {
            for (int j = 1; j < 3; j++) {
                double value = (double)scalars[indices[j]] - scalars[indices[0]];
                result = geomlib_vec3d_add(result, geomlib_vec3d_mul(triangle.gradients[j], value));
            }
        } else {
            for (int j = 1; j < 3; j++) {
                geomlib_vec3d value = geomlib_vec3d_sub(to_double(vectors[indices[j]]),
                                                      to_double(vectors[indices[0]]));
                geomlib_vec3d basis = triangle.gradients[j];
                if (operation == SURFACE_DIVERGENCE) {
                    result.x += geomlib_vec3d_dot(basis, value);
                } else {
                    result = geomlib_vec3d_add(result, geomlib_vec3d_cross(basis, value));
                }
            }
            if (operation == SURFACE_CURL) {
                result = geomlib_vec3d_mul(triangle.normal,
                    geomlib_vec3d_dot(result, triangle.normal));
            }
        }

        double weight = operation == SURFACE_LAPLACIAN ? triangle.area / 3.0 : triangle.area;
        for (int j = 0; j < 3; j++) {
            surface_accumulator *sum = &sums[indices[j]];
            sum->area += weight;
            if (operation == SURFACE_LAPLACIAN) {
                sum->value.x -= triangle.area * geomlib_vec3d_dot(triangle.gradients[j], result);
            } else {
                sum->value = geomlib_vec3d_add(sum->value, geomlib_vec3d_mul(result, weight));
            }
        }
    }

    for (size_t i = 0; i < count; i++) {
        geomlib_vec3d value = sums[i].area > 0 ? geomlib_vec3d_mul(sums[i].value, 1.0 / sums[i].area)
                                               : (geomlib_vec3d){0};
        if (!isfinite(value.x) || fabs(value.x) > FLT_MAX ||
            !isfinite(value.y) || fabs(value.y) > FLT_MAX ||
            !isfinite(value.z) || fabs(value.z) > FLT_MAX) {
            free(sums);
            return GEOMLIB_MESH_ERROR_INVALID;
        }
        sums[i].value = value;
    }
    for (size_t i = 0; i < count; i++) {
        if (scalar_out) scalar_out[i] = (float)sums[i].value.x;
        else vector_out[i] = to_float(sums[i].value);
    }
    free(sums);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_compute_gradient_scalar(const geomlib_trimesh *mesh,
                                    const geomlib_scalar_field *field,
                                    geomlib_gradient_method method,
                                    geomlib_vector_field *gradient) {
    if (!gradient || !valid_gradient_method(method) || !field ||
        !valid_surface(mesh, field->count) || !valid_scalar_values(field))
        return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_vector_field output = geomlib_vector_field_create(field->count);
    if (field->count && !output.values) return GEOMLIB_MESH_ERROR_MEMORY;
    int result = compute_surface_operator(mesh, field->count, field->values, NULL,
                                          SURFACE_GRADIENT, NULL, output.values);
    if (result == GEOMLIB_MESH_SUCCESS) *gradient = output;
    else geomlib_vector_field_free(&output);
    return result;
}

int geomlib_compute_divergence(const geomlib_trimesh *mesh, const geomlib_vector_field *field,
                               geomlib_scalar_field *divergence) {
    if (!field || !divergence || divergence->count != field->count ||
        (field->count && !divergence->values) || !valid_surface(mesh, field->count) ||
        !valid_vector_values(field)) return GEOMLIB_MESH_ERROR_INVALID;
    return compute_surface_operator(mesh, field->count, NULL, field->values,
                                    SURFACE_DIVERGENCE, divergence->values, NULL);
}

int geomlib_compute_curl(const geomlib_trimesh *mesh, const geomlib_vector_field *field,
                         geomlib_vector_field *curl) {
    if (!field || !curl || curl->count != field->count || (field->count && !curl->values) ||
        !valid_surface(mesh, field->count) || !valid_vector_values(field))
        return GEOMLIB_MESH_ERROR_INVALID;
    return compute_surface_operator(mesh, field->count, NULL, field->values,
                                    SURFACE_CURL, NULL, curl->values);
}

int geomlib_compute_laplacian(const geomlib_trimesh *mesh, const geomlib_scalar_field *field,
                              geomlib_scalar_field *laplacian) {
    if (!field || !laplacian || laplacian->count != field->count ||
        (field->count && !laplacian->values) || !valid_surface(mesh, field->count) ||
        !valid_scalar_values(field)) return GEOMLIB_MESH_ERROR_INVALID;
    return compute_surface_operator(mesh, field->count, field->values, NULL,
                                    SURFACE_LAPLACIAN, laplacian->values, NULL);
}

/* ============================================================================
 * Interpolation
 * ============================================================================ */

float geomlib_interpolate_face_value(
    float value_left,
    float value_right,
    float distance_ratio,
    geomlib_interpolation_scheme scheme
) {
    switch (scheme) {
        case GEOMLIB_INTERP_LINEAR:
            /* Linear interpolation */
            return value_left * (1.0f - distance_ratio) + value_right * distance_ratio;

        case GEOMLIB_INTERP_UPWIND:
            /* First-order upwind (use upstream value) */
            return (distance_ratio < 0.5f) ? value_left : value_right;

        case GEOMLIB_INTERP_CENTRAL:
            /* Central differencing (average) */
            return (value_left + value_right) * 0.5f;

        default:
            return value_left * (1.0f - distance_ratio) + value_right * distance_ratio;
    }
}

static bool interpolation_weights(const geomlib_trimesh *mesh, size_t count, geomlib_vec3 point,
                                  uint32_t indices[3], double weights[3]) {
    if (!finite_vector(point) || !valid_surface(mesh, count)) return false;

    double nearest = DBL_MAX;
    for (size_t face = 0; face < mesh->triangle_count; face++) {
        const uint32_t *corners = &mesh->indices[face * 3];
        geomlib_vec3 a = mesh->vertices[corners[0]], b = mesh->vertices[corners[1]],
                     c = mesh->vertices[corners[2]];
        surface_triangle triangle;
        if (!triangle_basis(a, b, c, &triangle)) continue;

        geomlib_vec3d origin = to_double(a);
        geomlib_vec3d edges[2] = {geomlib_vec3d_sub(to_double(b), origin),
                                  geomlib_vec3d_sub(to_double(c), origin)};
        geomlib_vec3d query = geomlib_vec3d_sub(to_double(point), origin);
        double face_weights[3];
        face_weights[1] = geomlib_vec3d_dot(triangle.gradients[1], query);
        face_weights[2] = geomlib_vec3d_dot(triangle.gradients[2], query);
        face_weights[0] = 1 - face_weights[1] - face_weights[2];
        if (face_weights[0] >= 0 && face_weights[1] >= 0 && face_weights[2] >= 0) {
            double distance = geomlib_vec3d_dot(triangle.normal, query);
            if (distance * distance < nearest) {
                nearest = distance * distance;
                memcpy(indices, corners, 3 * sizeof(uint32_t));
                memcpy(weights, face_weights, 3 * sizeof(double));
            }
            continue;
        }

        geomlib_vec3d vertices[3] = {{0}, edges[0], edges[1]};
        for (int i = 0; i < 3; i++) {
            int j = (i + 1) % 3;
            geomlib_vec3d edge = geomlib_vec3d_sub(vertices[j], vertices[i]);
            geomlib_vec3d offset = geomlib_vec3d_sub(query, vertices[i]);
            double t = geomlib_vec3d_dot(edge, offset) / geomlib_vec3d_dot(edge, edge);
            t = fmax(0, fmin(1, t));
            geomlib_vec3d delta = geomlib_vec3d_sub(offset, geomlib_vec3d_mul(edge, t));
            double distance = geomlib_vec3d_dot(delta, delta);
            if (distance < nearest) {
                nearest = distance;
                memcpy(indices, corners, 3 * sizeof(uint32_t));
                weights[i] = 1 - t;
                weights[j] = t;
                weights[(j + 1) % 3] = 0;
            }
        }
    }
    return nearest < DBL_MAX;
}

float geomlib_interpolate_scalar(const geomlib_trimesh *mesh, const geomlib_scalar_field *field,
                                 geomlib_vec3 point, geomlib_interpolation_scheme scheme) {
    (void)scheme;
    uint32_t indices[3];
    double weights[3];
    if (!field || !field->values || !interpolation_weights(mesh, field->count, point, indices, weights))
        return 0;
    double result = 0;
    for (int i = 0; i < 3; i++) result += weights[i] * field->values[indices[i]];
    return (float)result;
}

geomlib_vec3 geomlib_interpolate_vector(const geomlib_trimesh *mesh, const geomlib_vector_field *field,
                                        geomlib_vec3 point, geomlib_interpolation_scheme scheme) {
    (void)scheme;
    uint32_t indices[3];
    double weights[3];
    geomlib_vec3d result = {0};
    if (field && field->values && interpolation_weights(mesh, field->count, point, indices, weights)) {
        for (int i = 0; i < 3; i++) result = geomlib_vec3d_add(result,
            geomlib_vec3d_mul(to_double(field->values[indices[i]]), weights[i]));
    }
    return to_float(result);
}

/* ============================================================================
 * Streamline Tracing
 * ============================================================================ */

geomlib_streamline geomlib_streamline_create(size_t capacity) {
    geomlib_streamline stream;
    stream.points = malloc(capacity * sizeof(geomlib_vec3));
    stream.point_count = 0;
    stream.point_capacity = capacity;
    return stream;
}

void geomlib_streamline_free(geomlib_streamline *streamline) {
    if (!streamline) return;
    free(streamline->points);
    memset(streamline, 0, sizeof(geomlib_streamline));
}

/* Barycentric interpolation inside the containing triangle, with a nearest-vertex fallback. */
static geomlib_vec3 interpolate_velocity(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    geomlib_vec3 point
) {
    size_t limit = mesh->vertex_count < velocity->count ? mesh->vertex_count : velocity->count;
    for (size_t t = 0; mesh->indices && t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3], i1 = mesh->indices[t * 3 + 1],
                 i2 = mesh->indices[t * 3 + 2];
        if (i0 >= limit || i1 >= limit || i2 >= limit) continue;
        geomlib_vec3d p0 = {mesh->vertices[i0].x, mesh->vertices[i0].y, mesh->vertices[i0].z};
        geomlib_vec3d e1 = geomlib_vec3d_sub(
            (geomlib_vec3d){mesh->vertices[i1].x, mesh->vertices[i1].y, mesh->vertices[i1].z}, p0);
        geomlib_vec3d e2 = geomlib_vec3d_sub(
            (geomlib_vec3d){mesh->vertices[i2].x, mesh->vertices[i2].y, mesh->vertices[i2].z}, p0);
        geomlib_vec3d w = geomlib_vec3d_sub((geomlib_vec3d){point.x, point.y, point.z}, p0);
        double d00 = geomlib_vec3d_dot(e1, e1);
        double d01 = geomlib_vec3d_dot(e1, e2);
        double d11 = geomlib_vec3d_dot(e2, e2);
        double denominator = d00 * d11 - d01 * d01;
        if (!(fabs(denominator) > 1e-20)) continue;
        double d20 = geomlib_vec3d_dot(w, e1);
        double d21 = geomlib_vec3d_dot(w, e2);
        double v = (d11 * d20 - d01 * d21) / denominator;
        double wb = (d00 * d21 - d01 * d20) / denominator;
        double u = 1.0 - v - wb;
        const double tolerance = -1e-4;
        if (u < tolerance || v < tolerance || wb < tolerance) continue;
        geomlib_vec3d result = geomlib_vec3d_mul((geomlib_vec3d){
            velocity->values[i0].x, velocity->values[i0].y, velocity->values[i0].z}, u);
        result = geomlib_vec3d_add(result, geomlib_vec3d_mul((geomlib_vec3d){
            velocity->values[i1].x, velocity->values[i1].y, velocity->values[i1].z}, v));
        result = geomlib_vec3d_add(result, geomlib_vec3d_mul((geomlib_vec3d){
            velocity->values[i2].x, velocity->values[i2].y, velocity->values[i2].z}, wb));
        return (geomlib_vec3){(float)result.x, (float)result.y, (float)result.z};
    }

    float min_distance = FLT_MAX;
    geomlib_vec3 nearest = {0, 0, 0};
    for (size_t v = 0; v < limit; v++) {
        float distance = geomlib_vec3_length(geomlib_vec3_sub(point, mesh->vertices[v]));
        if (distance < min_distance) {
            min_distance = distance;
            nearest = velocity->values[v];
        }
    }
    return nearest;
}

int geomlib_trace_streamline(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    geomlib_vec3 seed,
    const geomlib_streamline_params *params,
    geomlib_streamline *streamline
) {
    if (!mesh || !velocity || !params || !streamline) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    streamline->point_count = 0;
    geomlib_vec3 current = seed;
    float total_length = 0.0f;

    for (int step = 0; step < params->max_steps; step++) {
        /* Add current point */
        if (streamline->point_count >= streamline->point_capacity) {
            streamline->point_capacity *= 2;
            streamline->points = realloc(streamline->points,
                streamline->point_capacity * sizeof(geomlib_vec3));
        }
        streamline->points[streamline->point_count++] = current;

        /* Check bounds */
        if (current.x < params->bounds.min.x || current.x > params->bounds.max.x ||
            current.y < params->bounds.min.y || current.y > params->bounds.max.y ||
            current.z < params->bounds.min.z || current.z > params->bounds.max.z) {
            break;  /* Exited domain */
        }

        /* Interpolate velocity */
        geomlib_vec3 vel = interpolate_velocity(mesh, velocity, current);

        /* Check for zero velocity */
        float vel_mag = geomlib_vec3_length(vel);
        if (vel_mag < 1e-6f) break;

        /* Integration step */
        geomlib_vec3 next;
        switch (params->method) {
            case GEOMLIB_INTEGRATE_EULER:
                /* Forward Euler */
                next = geomlib_vec3_add(current,
                    geomlib_vec3_mul(vel, params->step_size));
                break;

            case GEOMLIB_INTEGRATE_RK2: {
                /* RK2 / Midpoint method */
                geomlib_vec3 k1 = geomlib_vec3_mul(vel, params->step_size);
                geomlib_vec3 mid = geomlib_vec3_add(current,
                    geomlib_vec3_mul(k1, 0.5f));
                geomlib_vec3 vel_mid = interpolate_velocity(mesh, velocity, mid);
                geomlib_vec3 k2 = geomlib_vec3_mul(vel_mid, params->step_size);
                next = geomlib_vec3_add(current, k2);
                break;
            }

            case GEOMLIB_INTEGRATE_RK4: {
                /* RK4 */
                geomlib_vec3 k1 = geomlib_vec3_mul(vel, params->step_size);

                geomlib_vec3 p2 = geomlib_vec3_add(current, geomlib_vec3_mul(k1, 0.5f));
                geomlib_vec3 v2 = interpolate_velocity(mesh, velocity, p2);
                geomlib_vec3 k2 = geomlib_vec3_mul(v2, params->step_size);

                geomlib_vec3 p3 = geomlib_vec3_add(current, geomlib_vec3_mul(k2, 0.5f));
                geomlib_vec3 v3 = interpolate_velocity(mesh, velocity, p3);
                geomlib_vec3 k3 = geomlib_vec3_mul(v3, params->step_size);

                geomlib_vec3 p4 = geomlib_vec3_add(current, k3);
                geomlib_vec3 v4 = interpolate_velocity(mesh, velocity, p4);
                geomlib_vec3 k4 = geomlib_vec3_mul(v4, params->step_size);

                /* Weighted sum */
                geomlib_vec3 sum = geomlib_vec3_add(k1, geomlib_vec3_mul(k2, 2.0f));
                sum = geomlib_vec3_add(sum, geomlib_vec3_mul(k3, 2.0f));
                sum = geomlib_vec3_add(sum, k4);
                next = geomlib_vec3_add(current, geomlib_vec3_mul(sum, 1.0f / 6.0f));
                break;
            }

            default:
                next = geomlib_vec3_add(current, geomlib_vec3_mul(vel, params->step_size));
        }

        /* Update length */
        geomlib_vec3 step_vec = geomlib_vec3_sub(next, current);
        total_length += geomlib_vec3_length(step_vec);

        if (total_length > params->max_length) break;

        current = next;
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * CFD Quality Metrics
 * ============================================================================ */

geomlib_cfd_quality_thresholds geomlib_cfd_quality_thresholds_default(void) {
    geomlib_cfd_quality_thresholds thresholds;
    thresholds.max_skewness = 0.95f;
    thresholds.max_non_orthogonality = 70.0f;
    thresholds.max_expansion_ratio = 5.0f;
    return thresholds;
}

static geomlib_cfd_quality quality_from_fvm(const geomlib_fvm_mesh *fvm, size_t cell,
                                             const geomlib_cfd_quality_thresholds *thresholds) {
    geomlib_cfd_quality quality = {0};
    const geomlib_tetramesh *mesh = fvm->tetramesh;
    const uint32_t *vertices = mesh->tetrahedra[cell].vertices;
    double lengths[6], mean = 0.0;
    int edge = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            geomlib_vec3d a = to_double(mesh->vertices[vertices[i]]);
            geomlib_vec3d b = to_double(mesh->vertices[vertices[j]]);
            lengths[edge] = geomlib_vec3d_length(geomlib_vec3d_sub(b, a));
            mean += lengths[edge++];
        }
    }
    mean /= 6;
    double variance = 0;
    for (int i = 0; i < 6; i++) {
        double difference = lengths[i] - mean;
        variance += difference * difference;
    }
    quality.skewness = (float)fmin(1.0, sqrt(variance / 6) / mean);
    quality.expansion_ratio = 1;
    for (uint32_t i = 0; i < fvm->cell_face_count[cell]; i++) {
        uint32_t index = fvm->cell_faces[cell][i];
        const geomlib_fvm_face *face = &fvm->faces[index];
        if (face->cell_right == UINT32_MAX) continue;
        uint32_t other = face->cell_left == cell ? face->cell_right : face->cell_left;
        geomlib_vec3d delta = geomlib_vec3d_sub(to_double(fvm->cell_centers[other]),
                                              to_double(fvm->cell_centers[cell]));
        double distance = geomlib_vec3d_length(delta);
        double alignment = distance > 0 ? fabs(geomlib_vec3d_dot(to_double(face->normal), delta)) / distance : 0;
        float angle = (float)(acos(fmin(1, alignment)) * (180 / M_PI));
        quality.non_orthogonality = fmaxf(quality.non_orthogonality, angle);
        double volume = fvm->cell_volumes[cell], neighbor = fvm->cell_volumes[other];
        double ratio = fmax(volume, neighbor) / fmin(volume, neighbor);
        quality.expansion_ratio = fmaxf(quality.expansion_ratio, ratio > FLT_MAX ? INFINITY : (float)ratio);
    }
    quality.is_valid = quality.skewness <= thresholds->max_skewness &&
                       quality.non_orthogonality <= thresholds->max_non_orthogonality &&
                       quality.expansion_ratio <= thresholds->max_expansion_ratio;
    return quality;
}

geomlib_cfd_quality geomlib_compute_cfd_quality(const geomlib_tetramesh *mesh, uint32_t tet_idx,
                                               const geomlib_cfd_quality_thresholds *thresholds) {
    if (!mesh || !thresholds || tet_idx >= mesh->tetrahedron_count) return (geomlib_cfd_quality){0};
    geomlib_fvm_mesh fvm = {0};
    if (geomlib_fvm_build(&fvm, mesh) != GEOMLIB_MESH_SUCCESS) return (geomlib_cfd_quality){0};
    geomlib_cfd_quality quality = quality_from_fvm(&fvm, tet_idx, thresholds);
    geomlib_fvm_free(&fvm);
    return quality;
}

/* ============================================================================
 * Particle System
 * ============================================================================ */

geomlib_particle_system geomlib_particle_system_create(size_t capacity) {
    geomlib_particle_system system;
    system.particles = malloc(capacity * sizeof(geomlib_particle));
    system.particle_count = 0;
    system.particle_capacity = capacity;
    return system;
}

void geomlib_particle_system_free(geomlib_particle_system *system) {
    if (!system) return;
    free(system->particles);
    memset(system, 0, sizeof(geomlib_particle_system));
}

int geomlib_advect_particles(
    geomlib_particle_system *particles,
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    float dt,
    const geomlib_aabb *bounds
) {
    if (!particles || !mesh || !velocity) return GEOMLIB_MESH_ERROR_INVALID;

    for (size_t i = 0; i < particles->particle_count; i++) {
        if (!particles->particles[i].active) continue;

        geomlib_particle *p = &particles->particles[i];

        /* Interpolate velocity */
        geomlib_vec3 vel = interpolate_velocity(mesh, velocity, p->position);
        p->velocity = vel;

        /* Advect */
        geomlib_vec3 displacement = geomlib_vec3_mul(vel, dt);
        p->position = geomlib_vec3_add(p->position, displacement);
        p->time += dt;

        /* Check bounds */
        if (bounds) {
            if (p->position.x < bounds->min.x || p->position.x > bounds->max.x ||
                p->position.y < bounds->min.y || p->position.y > bounds->max.y ||
                p->position.z < bounds->min.z || p->position.z > bounds->max.z) {
                p->active = false;
            }
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

typedef struct cfd_tetrahedron {
    geomlib_vec3d corners[4];
    geomlib_vec3d cofactors[3];
    double determinant;
} cfd_tetrahedron;

static bool valid_volume_mesh(const geomlib_tetramesh *mesh) {
    return mesh && mesh->vertex_count <= UINT32_MAX &&
           mesh->tetrahedron_count <= UINT32_MAX / 4 &&
           (!mesh->vertex_count || mesh->vertices) &&
           (!mesh->tetrahedron_count || mesh->tetrahedra);
}

static int tetrahedron_basis(const geomlib_tetramesh *mesh, size_t cell, cfd_tetrahedron *tetra) {
    for (int i = 0; i < 4; i++) {
        uint32_t index = mesh->tetrahedra[cell].vertices[i];
        if (index >= mesh->vertex_count || !finite_vector(mesh->vertices[index]))
            return GEOMLIB_MESH_ERROR_INVALID;
        tetra->corners[i] = to_double(mesh->vertices[index]);
    }
    geomlib_vec3d edges[3];
    double scale = 1;
    for (int i = 0; i < 3; i++) {
        edges[i] = geomlib_vec3d_sub(tetra->corners[i + 1], tetra->corners[0]);
        scale *= geomlib_vec3d_length(edges[i]);
    }
    for (int i = 0; i < 3; i++)
        tetra->cofactors[i] = geomlib_vec3d_cross(edges[(i + 1) % 3], edges[(i + 2) % 3]);
    tetra->determinant = geomlib_vec3d_dot(edges[0], tetra->cofactors[0]);
    if (!(fabs(tetra->determinant) > DBL_EPSILON * scale)) return GEOMLIB_MESH_ERROR_DEGENERATE;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_compute_gradient_tetra(const geomlib_tetramesh *mesh, const geomlib_scalar_field *field,
                                   geomlib_gradient_method method, geomlib_vector_field *gradient) {
    if (!gradient || !valid_volume_mesh(mesh) || !valid_gradient_method(method) ||
        !field || field->count != mesh->vertex_count || !valid_scalar_values(field))
        return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_vector_field output = geomlib_vector_field_create(mesh->tetrahedron_count);
    if (output.count && !output.values) return GEOMLIB_MESH_ERROR_MEMORY;
    int result = GEOMLIB_MESH_SUCCESS;
    for (size_t cell = 0; cell < output.count; cell++) {
        cfd_tetrahedron tetra;
        result = tetrahedron_basis(mesh, cell, &tetra);
        if (result != GEOMLIB_MESH_SUCCESS) break;

        const uint32_t *indices = mesh->tetrahedra[cell].vertices;
        geomlib_vec3d value = {0};
        for (int i = 0; i < 3; i++) {
            double difference = (double)field->values[indices[i + 1]] - field->values[indices[0]];
            value = geomlib_vec3d_add(value, geomlib_vec3d_mul(tetra.cofactors[i], difference));
        }
        value = geomlib_vec3d_mul(value, 1.0 / tetra.determinant);
        if (fabs(value.x) > FLT_MAX || fabs(value.y) > FLT_MAX || fabs(value.z) > FLT_MAX) {
            result = GEOMLIB_MESH_ERROR_INVALID;
            break;
        }
        output.values[cell] = to_float(value);
    }
    if (result == GEOMLIB_MESH_SUCCESS) *gradient = output;
    else geomlib_vector_field_free(&output);
    return result;
}

typedef struct fvm_facet {
    uint32_t key[3];
    uint32_t cell;
    uint32_t opposite;
} fvm_facet;

static void sort_facet_key(uint32_t key[3]) {
    for (int i = 1; i < 3; i++) {
        for (int j = i; j > 0 && key[j - 1] > key[j]; j--) {
            uint32_t value = key[j - 1];
            key[j - 1] = key[j];
            key[j] = value;
        }
    }
}

static int compare_facets(const void *a, const void *b) {
    const fvm_facet *left = a, *right = b;
    for (int i = 0; i < 3; i++) {
        if (left->key[i] != right->key[i]) return left->key[i] < right->key[i] ? -1 : 1;
    }
    return (left->cell > right->cell) - (left->cell < right->cell);
}

static geomlib_vec3d facet_opposite_vertex(const geomlib_tetramesh *mesh, const fvm_facet *facet) {
    return to_double(mesh->vertices[mesh->tetrahedra[facet->cell].vertices[facet->opposite]]);
}

static int create_fvm_face(const geomlib_tetramesh *mesh, const fvm_facet *left,
                           const fvm_facet *right, geomlib_fvm_face *face) {
    memcpy(face->vertices, left->key, sizeof(left->key));
    geomlib_vec3d a = to_double(mesh->vertices[left->key[0]]);
    geomlib_vec3d b = to_double(mesh->vertices[left->key[1]]);
    geomlib_vec3d c = to_double(mesh->vertices[left->key[2]]);
    geomlib_vec3d normal = geomlib_vec3d_cross(geomlib_vec3d_sub(b, a), geomlib_vec3d_sub(c, a));
    double area2 = geomlib_vec3d_length(normal);
    if (!(area2 > 0) || area2 * 0.5 > FLT_MAX) return GEOMLIB_MESH_ERROR_DEGENERATE;
    if (geomlib_vec3d_dot(normal, geomlib_vec3d_sub(facet_opposite_vertex(mesh, left), a)) > 0) {
        face->vertices[1] = left->key[2];
        face->vertices[2] = left->key[1];
        normal = geomlib_vec3d_mul(normal, -1);
    }
    if (right && !(geomlib_vec3d_dot(normal,
                  geomlib_vec3d_sub(facet_opposite_vertex(mesh, right), a)) > 0))
        return GEOMLIB_MESH_ERROR_INVALID;

    face->vertex_count = 3;
    face->cell_left = left->cell;
    face->cell_right = right ? right->cell : UINT32_MAX;
    face->center = to_float(geomlib_vec3d_mul(geomlib_vec3d_add(geomlib_vec3d_add(a, b), c), 1.0 / 3));
    face->normal = to_float(geomlib_vec3d_mul(normal, 1.0 / area2));
    face->area = (float)(area2 * 0.5);
    return face->area > 0 ? GEOMLIB_MESH_SUCCESS : GEOMLIB_MESH_ERROR_DEGENERATE;
}

int geomlib_fvm_build(geomlib_fvm_mesh *fvm, const geomlib_tetramesh *mesh) {
    if (!fvm || !valid_volume_mesh(mesh)) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_fvm_mesh output = {0};
    output.tetramesh = (geomlib_tetramesh*)mesh;
    output.cell_count = mesh->tetrahedron_count;
    size_t cells = output.cell_count;
    if (!cells) {
        *fvm = output;
        return GEOMLIB_MESH_SUCCESS;
    }
    if (cells > SIZE_MAX / (4 * sizeof(fvm_facet)) ||
        cells > SIZE_MAX / (4 * sizeof(geomlib_fvm_face))) return GEOMLIB_MESH_ERROR_MEMORY;

    size_t facets = cells * 4;
    fvm_facet *entries = malloc(facets * sizeof(*entries));
    output.faces = calloc(facets, sizeof(*output.faces));
    output.cell_centers = malloc(cells * sizeof(*output.cell_centers));
    output.cell_volumes = malloc(cells * sizeof(*output.cell_volumes));
    output.cell_faces = calloc(cells, sizeof(*output.cell_faces));
    output.cell_face_count = malloc(cells * sizeof(*output.cell_face_count));
    if (output.cell_faces) output.cell_faces[0] = malloc(facets * sizeof(uint32_t));
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (!entries || !output.faces || !output.cell_centers || !output.cell_volumes ||
        !output.cell_faces || !output.cell_faces[0] || !output.cell_face_count) goto finish;

    for (size_t cell = 0; cell < cells; cell++) {
        cfd_tetrahedron tetra;
        result = tetrahedron_basis(mesh, cell, &tetra);
        if (result != GEOMLIB_MESH_SUCCESS) goto finish;
        double volume = fabs(tetra.determinant) / 6;
        if (volume > FLT_MAX || (float)volume == 0) {
            result = GEOMLIB_MESH_ERROR_DEGENERATE;
            goto finish;
        }
        output.cell_volumes[cell] = (float)volume;
        geomlib_vec3d center = {0};
        for (int j = 0; j < 4; j++) center = geomlib_vec3d_add(center, tetra.corners[j]);
        output.cell_centers[cell] = to_float(geomlib_vec3d_mul(center, 0.25));
        output.cell_faces[cell] = output.cell_faces[0] + cell * 4;
        output.cell_face_count[cell] = 4;
        for (uint32_t opposite = 0; opposite < 4; opposite++) {
            fvm_facet *entry = &entries[cell * 4 + opposite];
            entry->cell = (uint32_t)cell;
            entry->opposite = opposite;
            for (int k = 0; k < 3; k++)
                entry->key[k] = mesh->tetrahedra[cell].vertices[(opposite + k + 1) % 4];
            sort_facet_key(entry->key);
        }
    }

    qsort(entries, facets, sizeof(*entries), compare_facets);
    for (size_t start = 0; start < facets;) {
        size_t end = start + 1;
        while (end < facets && memcmp(entries[start].key, entries[end].key, sizeof(entries[0].key)) == 0)
            end++;
        if (end - start > 2) {
            result = GEOMLIB_MESH_ERROR_INVALID;
            goto finish;
        }
        const fvm_facet *right = end - start == 2 ? &entries[start + 1] : NULL;
        result = create_fvm_face(mesh, &entries[start], right, &output.faces[output.face_count]);
        if (result != GEOMLIB_MESH_SUCCESS) goto finish;
        for (size_t i = start; i < end; i++)
            output.cell_faces[entries[i].cell][entries[i].opposite] = (uint32_t)output.face_count;
        output.face_count++;
        start = end;
    }

finish:
    free(entries);
    if (result == GEOMLIB_MESH_SUCCESS) *fvm = output;
    else geomlib_fvm_free(&output);
    return result;
}

void geomlib_fvm_free(geomlib_fvm_mesh *fvm) {
    if (!fvm) return;
    free(fvm->faces);
    free(fvm->cell_centers);
    free(fvm->cell_volumes);
    if (fvm->cell_faces) free(fvm->cell_faces[0]);
    free(fvm->cell_faces);
    free(fvm->cell_face_count);
    memset(fvm, 0, sizeof(*fvm));
}

int geomlib_generate_streamlines(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    const geomlib_vec3 *seeds,
    size_t seed_count,
    const geomlib_streamline_params *params,
    geomlib_streamline *streamlines
) {
    if (!mesh || !velocity || !seeds || !streamlines) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    for (size_t i = 0; i < seed_count; i++) {
        geomlib_trace_streamline(mesh, velocity, seeds[i], params, &streamlines[i]);
    }

    return GEOMLIB_MESH_SUCCESS;
}

bool geomlib_check_cfd_mesh_quality(
    const geomlib_tetramesh *mesh,
    const geomlib_cfd_quality_thresholds *thresholds,
    size_t *failed_count
) {
    if (failed_count) *failed_count = mesh ? mesh->tetrahedron_count : 0;
    if (!mesh || !thresholds) return false;

    geomlib_fvm_mesh fvm = {0};
    if (geomlib_fvm_build(&fvm, mesh) != GEOMLIB_MESH_SUCCESS) return false;
    size_t failures = 0;
    for (size_t i = 0; i < fvm.cell_count; i++) {
        if (!quality_from_fvm(&fvm, i, thresholds).is_valid) failures++;
    }
    geomlib_fvm_free(&fvm);
    if (failed_count) *failed_count = failures;
    return failures == 0;
}
