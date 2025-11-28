/*
 * geomlib - Curvature Computation Implementation
 * Based on Meyer et al. 2003 discrete differential geometry operators
 */

#include "geomlib/geomlib_curvature.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Create curvature structure */
geomlib_mesh_curvature geomlib_curvature_create(size_t vertex_count) {
    geomlib_mesh_curvature curv;
    curv.curvatures = calloc(vertex_count, sizeof(geomlib_vertex_curvature));
    curv.vertex_areas = calloc(vertex_count, sizeof(float));
    curv.vertex_count = vertex_count;
    return curv;
}

/* Free curvature structure */
void geomlib_curvature_free(geomlib_mesh_curvature *curv) {
    if (!curv) return;
    if (curv->curvatures) {
        free(curv->curvatures);
        curv->curvatures = NULL;
    }
    if (curv->vertex_areas) {
        free(curv->vertex_areas);
        curv->vertex_areas = NULL;
    }
    curv->vertex_count = 0;
}

/* Compute cotangent weight for edge (v0, v1) with opposite vertex v2 */
float geomlib_cotan_weight(
    const geomlib_vec3 *v0,
    const geomlib_vec3 *v1,
    const geomlib_vec3 *v2
) {
    /* Compute vectors from v2 to v0 and v1 */
    geomlib_vec3 e0 = geomlib_vec3_sub(*v0, *v2);
    geomlib_vec3 e1 = geomlib_vec3_sub(*v1, *v2);

    /* cos(angle) = dot(e0, e1) / (|e0| * |e1|) */
    /* cot(angle) = cos(angle) / sin(angle) */
    /* sin(angle) = |cross(e0, e1)| / (|e0| * |e1|) */

    float dot = geomlib_vec3_dot(e0, e1);
    geomlib_vec3 cross = geomlib_vec3_cross(e0, e1);
    float cross_len = geomlib_vec3_length(cross);

    if (cross_len < 1e-10f) return 0.0f;  /* Degenerate triangle */

    return dot / cross_len;
}

/* Compute mixed Voronoi area for a vertex */
float geomlib_compute_voronoi_area(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
) {
    if (!mesh || vertex_idx >= mesh->vertex_count) return 0.0f;

    geomlib_vec3 v = mesh->vertices[vertex_idx];
    float area = 0.0f;

    /* Sum contributions from all adjacent triangles */
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        /* Check if this triangle contains our vertex */
        int vertex_pos = -1;
        if (i0 == vertex_idx) vertex_pos = 0;
        else if (i1 == vertex_idx) vertex_pos = 1;
        else if (i2 == vertex_idx) vertex_pos = 2;

        if (vertex_pos < 0) continue;

        /* Get triangle vertices */
        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        /* Compute triangle area */
        geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
        geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
        float tri_area = geomlib_vec3_length(geomlib_vec3_cross(e1, e2)) * 0.5f;

        /* Check if triangle is obtuse */
        geomlib_vec3 v0 = geomlib_vec3_sub(p1, p0);
        geomlib_vec3 v1 = geomlib_vec3_sub(p2, p0);
        geomlib_vec3 v2 = geomlib_vec3_sub(p2, p1);

        float dot0 = geomlib_vec3_dot(v0, v1);
        float dot1 = geomlib_vec3_dot(geomlib_vec3_negate(v0), v2);
        float dot2 = geomlib_vec3_dot(geomlib_vec3_negate(v1),
                                       geomlib_vec3_negate(v2));

        bool obtuse_at_vertex = (vertex_pos == 0 && dot0 < 0) ||
                                (vertex_pos == 1 && dot1 < 0) ||
                                (vertex_pos == 2 && dot2 < 0);

        bool obtuse = (dot0 < 0 || dot1 < 0 || dot2 < 0);

        if (obtuse) {
            /* Use special formula for obtuse triangles */
            if (obtuse_at_vertex) {
                area += tri_area * 0.5f;  /* Half the triangle area */
            } else {
                area += tri_area * 0.25f;  /* Quarter the triangle area */
            }
        } else {
            /* Use Voronoi area (1/8 sum of cot-weighted edge lengths) */
            /* For simplicity, use 1/3 of triangle area */
            area += tri_area / 3.0f;
        }
    }

    return area > 0.0f ? area : 1.0f;  /* Prevent division by zero */
}

/* Compute mean curvature using cotan-Laplacian */
int geomlib_compute_mean_curvature(
    const geomlib_trimesh *mesh,
    float *mean_curvatures
) {
    if (!mesh || !mean_curvatures) return GEOMLIB_MESH_ERROR_INVALID;

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Compute mean curvature for each vertex */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        geomlib_vec3 vertex = mesh->vertices[v];
        geomlib_vec3 laplacian = {0.0f, 0.0f, 0.0f};
        float weight_sum = 0.0f;

        /* Find all one-ring neighbors and compute cotan-weighted Laplacian */
        for (size_t he_idx = 0; he_idx < hemesh.halfedge_count; he_idx++) {
            if (hemesh.halfedges[he_idx].vertex != v) continue;

            uint32_t next = hemesh.halfedges[he_idx].next;
            if (next == UINT32_MAX) continue;

            uint32_t neighbor_idx = hemesh.halfedges[next].vertex;
            geomlib_vec3 neighbor = mesh->vertices[neighbor_idx];

            /* Find the two triangles sharing this edge */
            uint32_t face = hemesh.halfedges[he_idx].face;
            uint32_t twin = hemesh.halfedges[he_idx].twin;

            float cotan_sum = 0.0f;
            int cotan_count = 0;

            /* Cotan from current face */
            if (face != UINT32_MAX) {
                uint32_t i0 = mesh->indices[face * 3 + 0];
                uint32_t i1 = mesh->indices[face * 3 + 1];
                uint32_t i2 = mesh->indices[face * 3 + 2];

                /* Find the third vertex (opposite to edge) */
                uint32_t opposite = i0;
                if (opposite == v || opposite == neighbor_idx) opposite = i1;
                if (opposite == v || opposite == neighbor_idx) opposite = i2;

                geomlib_vec3 opp = mesh->vertices[opposite];
                float cotan = geomlib_cotan_weight(&vertex, &neighbor, &opp);
                cotan_sum += cotan;
                cotan_count++;
            }

            /* Cotan from twin face */
            if (twin != UINT32_MAX) {
                uint32_t twin_face = hemesh.halfedges[twin].face;
                if (twin_face != UINT32_MAX) {
                    uint32_t i0 = mesh->indices[twin_face * 3 + 0];
                    uint32_t i1 = mesh->indices[twin_face * 3 + 1];
                    uint32_t i2 = mesh->indices[twin_face * 3 + 2];

                    uint32_t opposite = i0;
                    if (opposite == v || opposite == neighbor_idx) opposite = i1;
                    if (opposite == v || opposite == neighbor_idx) opposite = i2;

                    geomlib_vec3 opp = mesh->vertices[opposite];
                    float cotan = geomlib_cotan_weight(&vertex, &neighbor, &opp);
                    cotan_sum += cotan;
                    cotan_count++;
                }
            }

            /* Add weighted edge to Laplacian */
            geomlib_vec3 edge = geomlib_vec3_sub(neighbor, vertex);
            geomlib_vec3 weighted = geomlib_vec3_mul(edge, cotan_sum);
            laplacian = geomlib_vec3_add(laplacian, weighted);
            weight_sum += cotan_sum;
        }

        /* Mean curvature H = |Laplacian| / (2 * area) */
        float voronoi_area = geomlib_compute_voronoi_area(mesh, v);
        float laplacian_mag = geomlib_vec3_length(laplacian);
        mean_curvatures[v] = laplacian_mag / (2.0f * voronoi_area);
    }

    geomlib_halfedge_free(&hemesh);
    return GEOMLIB_MESH_SUCCESS;
}

/* Compute Gaussian curvature using angle defect */
int geomlib_compute_gaussian_curvature(
    const geomlib_trimesh *mesh,
    float *gaussian_curvatures
) {
    if (!mesh || !gaussian_curvatures) return GEOMLIB_MESH_ERROR_INVALID;

    /* Initialize to 2*PI (for isolated vertices) */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        gaussian_curvatures[v] = 2.0f * M_PI;
    }

    /* Subtract angles from adjacent triangles */
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        /* Compute angles at each vertex */
        geomlib_vec3 e0 = geomlib_vec3_sub(p1, p0);
        geomlib_vec3 e1 = geomlib_vec3_sub(p2, p0);
        geomlib_vec3 e2 = geomlib_vec3_sub(p2, p1);

        float len_e0 = geomlib_vec3_length(e0);
        float len_e1 = geomlib_vec3_length(e1);
        float len_e2 = geomlib_vec3_length(e2);

        if (len_e0 < 1e-10f || len_e1 < 1e-10f || len_e2 < 1e-10f) continue;

        /* Angle at vertex 0 */
        float cos_angle0 = geomlib_vec3_dot(e0, e1) / (len_e0 * len_e1);
        cos_angle0 = fmaxf(-1.0f, fminf(1.0f, cos_angle0));
        float angle0 = acosf(cos_angle0);

        /* Angle at vertex 1 */
        geomlib_vec3 e0_neg = geomlib_vec3_negate(e0);
        float cos_angle1 = geomlib_vec3_dot(e0_neg, e2) / (len_e0 * len_e2);
        cos_angle1 = fmaxf(-1.0f, fminf(1.0f, cos_angle1));
        float angle1 = acosf(cos_angle1);

        /* Angle at vertex 2 */
        geomlib_vec3 e1_neg = geomlib_vec3_negate(e1);
        geomlib_vec3 e2_neg = geomlib_vec3_negate(e2);
        float cos_angle2 = geomlib_vec3_dot(e1_neg, e2_neg) / (len_e1 * len_e2);
        cos_angle2 = fmaxf(-1.0f, fminf(1.0f, cos_angle2));
        float angle2 = acosf(cos_angle2);

        /* Subtract from angle defect */
        gaussian_curvatures[i0] -= angle0;
        gaussian_curvatures[i1] -= angle1;
        gaussian_curvatures[i2] -= angle2;
    }

    /* Normalize by Voronoi area */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        float voronoi_area = geomlib_compute_voronoi_area(mesh, v);
        gaussian_curvatures[v] /= voronoi_area;
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* Compute principal curvatures from mean and Gaussian */
void geomlib_compute_principal_curvatures(
    float mean,
    float gaussian,
    float *k1,
    float *k2
) {
    /* k1, k2 are roots of: k^2 - 2*H*k + K = 0 */
    /* k1 = H + sqrt(H^2 - K) */
    /* k2 = H - sqrt(H^2 - K) */

    float discriminant = mean * mean - gaussian;

    if (discriminant < 0.0f) {
        /* Should not happen for valid surface, clamp to zero */
        discriminant = 0.0f;
    }

    float sqrt_disc = sqrtf(discriminant);
    *k1 = mean + sqrt_disc;
    *k2 = mean - sqrt_disc;
}

/* Compute all curvatures for mesh */
int geomlib_compute_curvature(
    const geomlib_trimesh *mesh,
    geomlib_mesh_curvature *curvature
) {
    if (!mesh || !curvature) return GEOMLIB_MESH_ERROR_INVALID;
    if (curvature->vertex_count != mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;

    /* Compute mean curvatures */
    float *mean_curv = malloc(mesh->vertex_count * sizeof(float));
    int result = geomlib_compute_mean_curvature(mesh, mean_curv);
    if (result != GEOMLIB_MESH_SUCCESS) {
        free(mean_curv);
        return result;
    }

    /* Compute Gaussian curvatures */
    float *gaussian_curv = malloc(mesh->vertex_count * sizeof(float));
    result = geomlib_compute_gaussian_curvature(mesh, gaussian_curv);
    if (result != GEOMLIB_MESH_SUCCESS) {
        free(mean_curv);
        free(gaussian_curv);
        return result;
    }

    /* Compute derived values for each vertex */
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        curvature->curvatures[v].mean = mean_curv[v];
        curvature->curvatures[v].gaussian = gaussian_curv[v];

        /* Compute principal curvatures */
        geomlib_compute_principal_curvatures(
            mean_curv[v],
            gaussian_curv[v],
            &curvature->curvatures[v].principal_max,
            &curvature->curvatures[v].principal_min
        );

        /* Store Voronoi area */
        curvature->vertex_areas[v] = geomlib_compute_voronoi_area(mesh, v);

        /* Compute principal directions */
        if (mesh->normals) {
            geomlib_vec3 normal = mesh->normals[v];

            /* Build tangent frame: find perpendicular vectors to normal */
            geomlib_vec3 tangent1, tangent2;

            /* Choose axis least aligned with normal */
            geomlib_vec3 axis;
            if (fabsf(normal.x) < fabsf(normal.y) && fabsf(normal.x) < fabsf(normal.z)) {
                axis = (geomlib_vec3){1.0f, 0.0f, 0.0f};
            } else if (fabsf(normal.y) < fabsf(normal.z)) {
                axis = (geomlib_vec3){0.0f, 1.0f, 0.0f};
            } else {
                axis = (geomlib_vec3){0.0f, 0.0f, 1.0f};
            }

            /* First tangent: perpendicular to normal */
            tangent1 = geomlib_vec3_normalize(geomlib_vec3_cross(normal, axis));

            /* Second tangent: complete orthonormal frame */
            tangent2 = geomlib_vec3_cross(normal, tangent1);

            /* Build 2x2 curvature tensor in tangent frame */
            float tensor[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
            float total_weight = 0.0f;

            /* Accumulate from incident edges */
            for (size_t t = 0; t < mesh->triangle_count; t++) {
                uint32_t i0 = mesh->indices[t * 3 + 0];
                uint32_t i1 = mesh->indices[t * 3 + 1];
                uint32_t i2 = mesh->indices[t * 3 + 2];

                /* Check if this triangle contains vertex v */
                if (i0 != v && i1 != v && i2 != v) continue;

                /* Get the other two vertices */
                uint32_t v1, v2;
                if (i0 == v) {
                    v1 = i1;
                    v2 = i2;
                } else if (i1 == v) {
                    v1 = i0;
                    v2 = i2;
                } else {
                    v1 = i0;
                    v2 = i1;
                }

                /* Edge vectors */
                geomlib_vec3 e1 = geomlib_vec3_sub(mesh->vertices[v1], mesh->vertices[v]);
                geomlib_vec3 e2 = geomlib_vec3_sub(mesh->vertices[v2], mesh->vertices[v]);

                /* Project edges onto tangent plane */
                geomlib_vec3 e1_proj = geomlib_vec3_sub(e1, geomlib_vec3_mul(normal, geomlib_vec3_dot(e1, normal)));
                geomlib_vec3 e2_proj = geomlib_vec3_sub(e2, geomlib_vec3_mul(normal, geomlib_vec3_dot(e2, normal)));

                /* Get coordinates in tangent frame */
                float e1_t1 = geomlib_vec3_dot(e1_proj, tangent1);
                float e1_t2 = geomlib_vec3_dot(e1_proj, tangent2);
                float e2_t1 = geomlib_vec3_dot(e2_proj, tangent1);
                float e2_t2 = geomlib_vec3_dot(e2_proj, tangent2);

                /* Normal curvature in edge direction */
                float len_e1 = geomlib_vec3_length(e1);
                float len_e2 = geomlib_vec3_length(e2);

                if (len_e1 > 1e-10f && len_e2 > 1e-10f) {
                    /* Weight by edge length */
                    float weight = len_e1 + len_e2;

                    /* Accumulate to tensor (simplified) */
                    tensor[0][0] += weight * e1_t1 * e1_t1;
                    tensor[0][1] += weight * e1_t1 * e1_t2;
                    tensor[1][0] += weight * e2_t1 * e2_t2;
                    tensor[1][1] += weight * e2_t2 * e2_t2;

                    total_weight += weight;
                }
            }

            if (total_weight > 1e-10f) {
                /* Normalize tensor */
                tensor[0][0] /= total_weight;
                tensor[0][1] /= total_weight;
                tensor[1][0] /= total_weight;
                tensor[1][1] /= total_weight;

                /* Compute eigenvalues and eigenvectors of 2x2 symmetric matrix */
                /* For symmetric 2x2: [[a, b], [b, c]] */
                float a = tensor[0][0];
                float b = (tensor[0][1] + tensor[1][0]) * 0.5f;  /* Symmetrize */
                float c = tensor[1][1];

                /* Eigenvalues: lambda = (a+c)/2 ± sqrt((a-c)^2/4 + b^2) */
                float trace = a + c;
                float det = a * c - b * b;
                float disc = trace * trace * 0.25f - det;

                if (disc >= 0.0f) {
                    float sqrt_disc = sqrtf(disc);
                    float lambda1 = trace * 0.5f + sqrt_disc;
                    float lambda2 = trace * 0.5f - sqrt_disc;

                    /* Eigenvector for lambda1: direction of maximum curvature */
                    geomlib_vec3 dir_max;
                    if (fabsf(b) > 1e-10f) {
                        float v1_x = lambda1 - c;
                        float v1_y = b;
                        float len = sqrtf(v1_x * v1_x + v1_y * v1_y);
                        v1_x /= len;
                        v1_y /= len;
                        dir_max = geomlib_vec3_add(
                            geomlib_vec3_mul(tangent1, v1_x),
                            geomlib_vec3_mul(tangent2, v1_y)
                        );
                    } else {
                        /* Diagonal matrix - principal directions align with tangent frame */
                        dir_max = (lambda1 > lambda2) ? tangent1 : tangent2;
                    }

                    /* Minimum direction is perpendicular in tangent plane */
                    geomlib_vec3 dir_min = geomlib_vec3_cross(normal, dir_max);

                    curvature->curvatures[v].principal_dir_max = geomlib_vec3_normalize(dir_max);
                    curvature->curvatures[v].principal_dir_min = geomlib_vec3_normalize(dir_min);
                } else {
                    /* Fallback to tangent frame */
                    curvature->curvatures[v].principal_dir_max = tangent1;
                    curvature->curvatures[v].principal_dir_min = tangent2;
                }
            } else {
                /* No incident triangles or degenerate - use tangent frame */
                curvature->curvatures[v].principal_dir_max = tangent1;
                curvature->curvatures[v].principal_dir_min = tangent2;
            }
        } else {
            /* No normals available - cannot compute directions */
            curvature->curvatures[v].principal_dir_max = (geomlib_vec3){0.0f, 0.0f, 0.0f};
            curvature->curvatures[v].principal_dir_min = (geomlib_vec3){0.0f, 0.0f, 0.0f};
        }
    }

    free(mean_curv);
    free(gaussian_curv);

    return GEOMLIB_MESH_SUCCESS;
}
