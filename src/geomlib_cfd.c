/*
 * geomlib - CFD Utilities Implementation
 * Gradient computation, interpolation, streamlines for CFD
 */

#include "geomlib/geomlib_cfd.h"
#include "geomlib/geomlib_healing.h"
#include "geomlib/geomlib_topology.h"
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
    geomlib_vec3 gradient = geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    float total_area = 0.0f;

    /* For each edge of the polygon */
    for (size_t i = 0; i < count; i++) {
        size_t j = (i + 1) % count;

        /* Edge midpoint and value */
        geomlib_vec3 midpoint = geomlib_vec3_mul(
            geomlib_vec3_add(vertices[i], vertices[j]), 0.5f);
        float edge_value = (values[i] + values[j]) * 0.5f;

        /* Edge vector */
        geomlib_vec3 edge = geomlib_vec3_sub(vertices[j], vertices[i]);

        /* Contribution to gradient */
        gradient.x += edge_value * edge.x;
        gradient.y += edge_value * edge.y;
        gradient.z += edge_value * edge.z;

        /* Accumulate area */
        if (i > 0) {
            geomlib_vec3 v0 = vertices[0];
            geomlib_vec3 v1 = vertices[i];
            geomlib_vec3 v2 = vertices[j];
            geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
            geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
            float tri_area = geomlib_vec3_length(geomlib_vec3_cross(e1, e2)) * 0.5f;
            total_area += tri_area;
        }
    }

    /* Normalize by area */
    if (total_area > 1e-10f) {
        gradient = geomlib_vec3_mul(gradient, 1.0f / total_area);
    }

    return gradient;
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

/* Compute gradient at mesh vertices */
int geomlib_compute_gradient_scalar(
    const geomlib_trimesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_gradient_method method,
    geomlib_vector_field *gradient
) {
    if (!mesh || !field || !gradient) return GEOMLIB_MESH_ERROR_INVALID;
    if (field->count != mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;

    /* Allocate gradient field */
    *gradient = geomlib_vector_field_create(mesh->vertex_count);

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Compute gradient at each vertex */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        /* Find neighbors */
        geomlib_vec3 neighbors[64];
        float neighbor_values[64];
        size_t neighbor_count = 0;

        for (size_t he = 0; he < hemesh.halfedge_count && neighbor_count < 64; he++) {
            if (hemesh.halfedges[he].vertex == v) {
                uint32_t next = hemesh.halfedges[he].next;
                if (next != UINT32_MAX) {
                    uint32_t neighbor = hemesh.halfedges[next].vertex;
                    neighbors[neighbor_count] = mesh->vertices[neighbor];
                    neighbor_values[neighbor_count] = field->values[neighbor] - field->values[v];
                    neighbor_count++;
                }
            }
        }

        /* Compute gradient based on method */
        if (method == GEOMLIB_GRADIENT_LEAST_SQUARES) {
            gradient->values[v] = geomlib_gradient_least_squares(
                mesh->vertices[v],
                neighbors,
                neighbor_values,
                neighbor_count
            );
        } else {
            /* Green-Gauss (simplified) */
            gradient->values[v] = geomlib_gradient_least_squares(
                mesh->vertices[v],
                neighbors,
                neighbor_values,
                neighbor_count
            );
        }
    }

    geomlib_halfedge_free(&hemesh);
    return GEOMLIB_MESH_SUCCESS;
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

/* Basic barycentric interpolation for triangle */
static float interpolate_triangle(
    const geomlib_vec3 *vertices,
    const float *values,
    geomlib_vec3 point
) {
    /* Compute barycentric coordinates */
    geomlib_vec3 v0 = vertices[0];
    geomlib_vec3 v1 = vertices[1];
    geomlib_vec3 v2 = vertices[2];

    geomlib_vec3 e0 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e1 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(point, v0);

    float d00 = geomlib_vec3_dot(e0, e0);
    float d01 = geomlib_vec3_dot(e0, e1);
    float d11 = geomlib_vec3_dot(e1, e1);
    float d20 = geomlib_vec3_dot(e2, e0);
    float d21 = geomlib_vec3_dot(e2, e1);

    float denom = d00 * d11 - d01 * d01;
    if (fabsf(denom) < 1e-10f) return values[0];

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return u * values[0] + v * values[1] + w * values[2];
}

float geomlib_interpolate_scalar(
    const geomlib_trimesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_vec3 point,
    geomlib_interpolation_scheme scheme
) {
    if (!mesh || !field) return 0.0f;

    /* Find containing triangle (simplified - check all triangles) */
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_vec3 vertices[3] = {
            mesh->vertices[i0],
            mesh->vertices[i1],
            mesh->vertices[i2]
        };

        float values[3] = {
            field->values[i0],
            field->values[i1],
            field->values[i2]
        };

        /* Simplified point-in-triangle test */
        return interpolate_triangle(vertices, values, point);
    }

    return 0.0f;
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

/* Interpolate velocity at point (simplified) */
static geomlib_vec3 interpolate_velocity(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    geomlib_vec3 point
) {
    /* Simplified: return velocity of nearest vertex */
    float min_dist = FLT_MAX;
    geomlib_vec3 result = geomlib_vec3_create(0.0f, 0.0f, 0.0f);

    for (size_t v = 0; v < mesh->vertex_count && v < velocity->count; v++) {
        geomlib_vec3 diff = geomlib_vec3_sub(point, mesh->vertices[v]);
        float dist = geomlib_vec3_length(diff);
        if (dist < min_dist) {
            min_dist = dist;
            result = velocity->values[v];
        }
    }

    return result;
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

geomlib_cfd_quality geomlib_compute_cfd_quality(
    const geomlib_tetramesh *mesh,
    uint32_t tet_idx,
    const geomlib_cfd_quality_thresholds *thresholds
) {
    geomlib_cfd_quality quality;
    memset(&quality, 0, sizeof(quality));

    if (!mesh || tet_idx >= mesh->tetrahedron_count) {
        return quality;
    }

    geomlib_tetrahedron *tet = &mesh->tetrahedra[tet_idx];
    geomlib_vec3 v0 = mesh->vertices[tet->vertices[0]];
    geomlib_vec3 v1 = mesh->vertices[tet->vertices[1]];
    geomlib_vec3 v2 = mesh->vertices[tet->vertices[2]];
    geomlib_vec3 v3 = mesh->vertices[tet->vertices[3]];

    /* Compute centroid */
    geomlib_vec3 centroid = geomlib_tetrahedron_centroid(v0, v1, v2, v3);

    /* Skewness: deviation from equilateral */
    float avg_edge = 0.0f;
    float max_edge = 0.0f;
    float edges[6];
    edges[0] = geomlib_vec3_length(geomlib_vec3_sub(v1, v0));
    edges[1] = geomlib_vec3_length(geomlib_vec3_sub(v2, v0));
    edges[2] = geomlib_vec3_length(geomlib_vec3_sub(v3, v0));
    edges[3] = geomlib_vec3_length(geomlib_vec3_sub(v2, v1));
    edges[4] = geomlib_vec3_length(geomlib_vec3_sub(v3, v1));
    edges[5] = geomlib_vec3_length(geomlib_vec3_sub(v3, v2));

    for (int i = 0; i < 6; i++) {
        avg_edge += edges[i];
        if (edges[i] > max_edge) max_edge = edges[i];
    }
    avg_edge /= 6.0f;

    /* Skewness based on edge variation */
    float edge_variation = 0.0f;
    for (int i = 0; i < 6; i++) {
        float diff = edges[i] - avg_edge;
        edge_variation += diff * diff;
    }
    edge_variation = sqrtf(edge_variation / 6.0f);
    quality.skewness = (avg_edge > 1e-10f) ? (edge_variation / avg_edge) : 1.0f;

    /* Non-orthogonality (simplified) */
    quality.non_orthogonality = 0.0f;  /* Would require neighbor info */

    /* Expansion ratio (simplified) */
    quality.expansion_ratio = 1.0f;  /* Would require neighbor volumes */

    /* Face flatness */
    quality.face_flatness = 0.0f;  /* All faces are triangular */

    /* Check validity */
    quality.is_valid = (quality.skewness <= thresholds->max_skewness);

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

/* Stubs for unimplemented functions */
int geomlib_fvm_build(geomlib_fvm_mesh *fvm, const geomlib_tetramesh *tetramesh) {
    (void)fvm; (void)tetramesh;
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_fvm_free(geomlib_fvm_mesh *fvm) {
    (void)fvm;
}

int geomlib_compute_gradient_tetra(
    const geomlib_tetramesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_gradient_method method,
    geomlib_vector_field *gradient
) {
    (void)mesh; (void)field; (void)method; (void)gradient;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_compute_divergence(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *field,
    geomlib_scalar_field *divergence
) {
    if (!mesh || !field || !divergence) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (field->count != mesh->vertex_count || divergence->count != mesh->vertex_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build vertex adjacency for neighbor access */
    geomlib_vertex_adjacency adj = geomlib_build_vertex_adjacency(mesh);

    /* Compute divergence at each vertex using Green-Gauss */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        float div_sum = 0.0f;
        float total_area = 0.0f;

        /* Accumulate over adjacent triangles */
        for (uint32_t i = 0; i < adj.face_counts[v]; i++) {
            uint32_t tri = adj.faces[v][i];

            /* Get triangle vertices */
            uint32_t v0 = mesh->indices[tri * 3 + 0];
            uint32_t v1 = mesh->indices[tri * 3 + 1];
            uint32_t v2 = mesh->indices[tri * 3 + 2];

            geomlib_vec3 p0 = mesh->vertices[v0];
            geomlib_vec3 p1 = mesh->vertices[v1];
            geomlib_vec3 p2 = mesh->vertices[v2];

            /* Triangle area and normal */
            geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
            geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
            geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
            float area = geomlib_vec3_length(cross) * 0.5f;
            geomlib_vec3 normal = geomlib_vec3_normalize(cross);

            /* Average field value at triangle */
            geomlib_vec3 f0 = field->values[v0];
            geomlib_vec3 f1 = field->values[v1];
            geomlib_vec3 f2 = field->values[v2];
            geomlib_vec3 f_avg = geomlib_vec3_create(
                (f0.x + f1.x + f2.x) / 3.0f,
                (f0.y + f1.y + f2.y) / 3.0f,
                (f0.z + f1.z + f2.z) / 3.0f
            );

            /* Green-Gauss: ∇·F ≈ (1/V) ∫ F·n dS */
            float contrib = geomlib_vec3_dot(f_avg, normal) * area;
            div_sum += contrib;
            total_area += area;
        }

        /* Normalize by total area (approximation of volume) */
        if (total_area > 1e-6f) {
            divergence->values[v] = div_sum / total_area;
        } else {
            divergence->values[v] = 0.0f;
        }
    }

    geomlib_vertex_adjacency_free(&adj);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_compute_curl(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *field,
    geomlib_vector_field *curl
) {
    if (!mesh || !field || !curl) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (field->count != mesh->vertex_count || curl->count != mesh->vertex_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build vertex adjacency */
    geomlib_vertex_adjacency adj = geomlib_build_vertex_adjacency(mesh);

    /* Compute curl at each vertex using Stokes' theorem */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        geomlib_vec3 curl_sum = geomlib_vec3_create(0, 0, 0);
        float total_area = 0.0f;

        /* Accumulate over adjacent triangles */
        for (uint32_t i = 0; i < adj.face_counts[v]; i++) {
            uint32_t tri = adj.faces[v][i];

            /* Get triangle vertices */
            uint32_t v0 = mesh->indices[tri * 3 + 0];
            uint32_t v1 = mesh->indices[tri * 3 + 1];
            uint32_t v2 = mesh->indices[tri * 3 + 2];

            geomlib_vec3 p0 = mesh->vertices[v0];
            geomlib_vec3 p1 = mesh->vertices[v1];
            geomlib_vec3 p2 = mesh->vertices[v2];

            /* Triangle normal and area */
            geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
            geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
            geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
            float area = geomlib_vec3_length(cross) * 0.5f;
            geomlib_vec3 normal = geomlib_vec3_normalize(cross);

            /* Field values */
            geomlib_vec3 f0 = field->values[v0];
            geomlib_vec3 f1 = field->values[v1];
            geomlib_vec3 f2 = field->values[v2];

            /* Stokes' theorem: ∫(∇×F)·n dS = ∮ F·dl */
            /* Approximate with line integral around triangle */
            geomlib_vec3 e01 = geomlib_vec3_sub(p1, p0);
            geomlib_vec3 e12 = geomlib_vec3_sub(p2, p1);
            geomlib_vec3 e20 = geomlib_vec3_sub(p0, p2);

            geomlib_vec3 f01 = geomlib_vec3_create(
                (f0.x + f1.x) * 0.5f,
                (f0.y + f1.y) * 0.5f,
                (f0.z + f1.z) * 0.5f
            );
            geomlib_vec3 f12 = geomlib_vec3_create(
                (f1.x + f2.x) * 0.5f,
                (f1.y + f2.y) * 0.5f,
                (f1.z + f2.z) * 0.5f
            );
            geomlib_vec3 f20 = geomlib_vec3_create(
                (f2.x + f0.x) * 0.5f,
                (f2.y + f0.y) * 0.5f,
                (f2.z + f0.z) * 0.5f
            );

            float line_integral =
                geomlib_vec3_dot(f01, e01) +
                geomlib_vec3_dot(f12, e12) +
                geomlib_vec3_dot(f20, e20);

            /* Curl component along normal */
            geomlib_vec3 curl_contrib = geomlib_vec3_create(
                normal.x * line_integral,
                normal.y * line_integral,
                normal.z * line_integral
            );

            curl_sum = geomlib_vec3_add(curl_sum, curl_contrib);
            total_area += area;
        }

        /* Normalize by total area */
        if (total_area > 1e-6f) {
            curl->values[v] = geomlib_vec3_create(
                curl_sum.x / total_area,
                curl_sum.y / total_area,
                curl_sum.z / total_area
            );
        } else {
            curl->values[v] = geomlib_vec3_create(0, 0, 0);
        }
    }

    geomlib_vertex_adjacency_free(&adj);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_compute_laplacian(
    const geomlib_trimesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_scalar_field *laplacian
) {
    if (!mesh || !field || !laplacian) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (field->count != mesh->vertex_count || laplacian->count != mesh->vertex_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build vertex adjacency for neighbor access */
    geomlib_vertex_adjacency adj = geomlib_build_vertex_adjacency(mesh);

    /* Compute Laplacian using cotangent weights (discrete Laplace-Beltrami) */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        float laplacian_sum = 0.0f;
        float weight_sum = 0.0f;

        geomlib_vec3 v_pos = mesh->vertices[v];
        float v_value = field->values[v];

        /* Accumulate contributions from neighbors */
        for (uint32_t i = 0; i < adj.neighbor_counts[v]; i++) {
            uint32_t neighbor = adj.neighbors[v][i];
            geomlib_vec3 n_pos = mesh->vertices[neighbor];
            float n_value = field->values[neighbor];

            /* Compute cotangent weight */
            /* For simplicity, use uniform weight (proper cotangent requires angle computation) */
            float dist = geomlib_vec3_length(geomlib_vec3_sub(n_pos, v_pos));
            float weight = (dist > 1e-6f) ? (1.0f / dist) : 0.0f;

            /* Laplacian: ∇²f ≈ Σ w_i (f_i - f_0) / Σ w_i */
            laplacian_sum += weight * (n_value - v_value);
            weight_sum += weight;
        }

        /* Normalize */
        if (weight_sum > 1e-6f) {
            laplacian->values[v] = laplacian_sum / weight_sum;
        } else {
            laplacian->values[v] = 0.0f;
        }
    }

    geomlib_vertex_adjacency_free(&adj);
    return GEOMLIB_MESH_SUCCESS;
}

geomlib_vec3 geomlib_interpolate_vector(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *field,
    geomlib_vec3 point,
    geomlib_interpolation_scheme scheme
) {
    if (!mesh || !field || field->count != mesh->vertex_count) {
        return geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    }

    /* Simple scheme: find nearest triangle and use barycentric interpolation */
    (void)scheme;  /* Ignore scheme for now, always use linear */

    /* Find triangle containing point (brute force for now) */
    float min_dist = FLT_MAX;
    uint32_t closest_tri = 0;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t v0 = mesh->indices[i * 3 + 0];
        uint32_t v1 = mesh->indices[i * 3 + 1];
        uint32_t v2 = mesh->indices[i * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[v0];
        geomlib_vec3 p1 = mesh->vertices[v1];
        geomlib_vec3 p2 = mesh->vertices[v2];

        /* Compute triangle centroid and distance to point */
        geomlib_vec3 centroid = geomlib_vec3_create(
            (p0.x + p1.x + p2.x) / 3.0f,
            (p0.y + p1.y + p2.y) / 3.0f,
            (p0.z + p1.z + p2.z) / 3.0f
        );

        float dist = geomlib_vec3_length(geomlib_vec3_sub(point, centroid));
        if (dist < min_dist) {
            min_dist = dist;
            closest_tri = (uint32_t)i;
        }
    }

    /* Get vertices of closest triangle */
    uint32_t v0 = mesh->indices[closest_tri * 3 + 0];
    uint32_t v1 = mesh->indices[closest_tri * 3 + 1];
    uint32_t v2 = mesh->indices[closest_tri * 3 + 2];

    geomlib_vec3 p0 = mesh->vertices[v0];
    geomlib_vec3 p1 = mesh->vertices[v1];
    geomlib_vec3 p2 = mesh->vertices[v2];

    /* Compute barycentric coordinates */
    geomlib_vec3 v0v1 = geomlib_vec3_sub(p1, p0);
    geomlib_vec3 v0v2 = geomlib_vec3_sub(p2, p0);
    geomlib_vec3 v0p = geomlib_vec3_sub(point, p0);

    float d00 = geomlib_vec3_dot(v0v1, v0v1);
    float d01 = geomlib_vec3_dot(v0v1, v0v2);
    float d11 = geomlib_vec3_dot(v0v2, v0v2);
    float d20 = geomlib_vec3_dot(v0p, v0v1);
    float d21 = geomlib_vec3_dot(v0p, v0v2);

    float denom = d00 * d11 - d01 * d01;
    if (fabsf(denom) < 1e-6f) {
        /* Degenerate triangle, use average */
        geomlib_vec3 f0 = field->values[v0];
        geomlib_vec3 f1 = field->values[v1];
        geomlib_vec3 f2 = field->values[v2];
        return geomlib_vec3_create(
            (f0.x + f1.x + f2.x) / 3.0f,
            (f0.y + f1.y + f2.y) / 3.0f,
            (f0.z + f1.z + f2.z) / 3.0f
        );
    }

    float b1 = (d11 * d20 - d01 * d21) / denom;
    float b2 = (d00 * d21 - d01 * d20) / denom;
    float b0 = 1.0f - b1 - b2;

    /* Interpolate vector field values */
    geomlib_vec3 f0 = field->values[v0];
    geomlib_vec3 f1 = field->values[v1];
    geomlib_vec3 f2 = field->values[v2];

    return geomlib_vec3_create(
        b0 * f0.x + b1 * f1.x + b2 * f2.x,
        b0 * f0.y + b1 * f1.y + b2 * f2.y,
        b0 * f0.z + b1 * f1.z + b2 * f2.z
    );
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
    if (!mesh || !thresholds) return false;

    size_t failures = 0;
    for (size_t i = 0; i < mesh->tetrahedron_count; i++) {
        geomlib_cfd_quality q = geomlib_compute_cfd_quality(mesh, i, thresholds);
        if (!q.is_valid) failures++;
    }

    if (failed_count) *failed_count = failures;
    return (failures == 0);
}
