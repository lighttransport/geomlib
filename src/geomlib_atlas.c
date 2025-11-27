/*
 * geomlib - Texture Atlas Generation Implementation
 */

#include "geomlib/geomlib_atlas.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Default parameters */
geomlib_atlas_params geomlib_atlas_params_default(void) {
    geomlib_atlas_params params;
    params.segmentation = GEOMLIB_ATLAS_SEG_AUTO;
    params.packing = GEOMLIB_ATLAS_PACK_SHELF;
    params.parameterization = GEOMLIB_PARAM_LSCM;

    params.max_chart_area = 0.2f;         /* Max 20% of mesh per chart */
    params.normal_threshold = 60.0f;      /* 60 degree threshold */
    params.distortion_threshold = 2.0f;   /* Max 2x stretch */

    params.gutter_width = 0.005f;         /* 0.5% gutter */
    params.atlas_resolution = 1024;
    params.normalize_charts = true;

    params.max_iterations = 10;
    params.optimize_rotation = true;

    return params;
}

/* Create atlas */
geomlib_atlas geomlib_atlas_create(size_t max_charts) {
    geomlib_atlas atlas;
    atlas.charts = malloc(max_charts * sizeof(geomlib_atlas_chart));
    atlas.chart_count = 0;
    atlas.uvs = NULL;
    atlas.uv_count = 0;
    atlas.coverage = 0.0f;
    atlas.avg_stretch = 1.0f;
    atlas.max_stretch = 1.0f;
    return atlas;
}

/* Free atlas */
void geomlib_atlas_free(geomlib_atlas *atlas) {
    if (!atlas) return;

    for (size_t i = 0; i < atlas->chart_count; i++) {
        free(atlas->charts[i].triangles);
        free(atlas->charts[i].uvs);
        free(atlas->charts[i].vertex_map);
    }
    free(atlas->charts);
    free(atlas->uvs);

    memset(atlas, 0, sizeof(geomlib_atlas));
}

/* Helper: Triangle normal */
static geomlib_vec3 compute_triangle_normal(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 n = geomlib_vec3_cross(e1, e2);
    return geomlib_vec3_normalize(n);
}

/* Helper: Triangle area */
static float compute_triangle_area(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
    return geomlib_vec3_length(cross) * 0.5f;
}

/* Segment mesh by normal clustering */
static int segment_by_normals(
    const geomlib_trimesh *mesh,
    const geomlib_atlas_params *params,
    geomlib_atlas *atlas
) {
    if (!mesh || !atlas) return GEOMLIB_MESH_ERROR_INVALID;

    /* Compute triangle normals */
    geomlib_vec3 *tri_normals = malloc(mesh->triangle_count * sizeof(geomlib_vec3));
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];
        tri_normals[t] = compute_triangle_normal(
            mesh->vertices[i0],
            mesh->vertices[i1],
            mesh->vertices[i2]
        );
    }

    /* Simple greedy clustering */
    bool *assigned = calloc(mesh->triangle_count, sizeof(bool));
    float cos_threshold = cosf(params->normal_threshold * M_PI / 180.0f);

    atlas->chart_count = 0;

    for (size_t seed = 0; seed < mesh->triangle_count; seed++) {
        if (assigned[seed]) continue;

        /* Start new chart */
        geomlib_atlas_chart *chart = &atlas->charts[atlas->chart_count];
        chart->triangles = malloc(mesh->triangle_count * sizeof(uint32_t));
        chart->triangle_count = 0;

        /* Grow chart from seed */
        uint32_t *queue = malloc(mesh->triangle_count * sizeof(uint32_t));
        size_t queue_size = 0;
        queue[queue_size++] = seed;
        assigned[seed] = true;

        geomlib_vec3 chart_normal = tri_normals[seed];

        while (queue_size > 0) {
            uint32_t current = queue[--queue_size];
            chart->triangles[chart->triangle_count++] = current;

            /* Find adjacent triangles with similar normals */
            uint32_t i0 = mesh->indices[current * 3 + 0];
            uint32_t i1 = mesh->indices[current * 3 + 1];
            uint32_t i2 = mesh->indices[current * 3 + 2];

            for (size_t t = 0; t < mesh->triangle_count; t++) {
                if (assigned[t]) continue;

                /* Check if triangles share an edge */
                uint32_t ti0 = mesh->indices[t * 3 + 0];
                uint32_t ti1 = mesh->indices[t * 3 + 1];
                uint32_t ti2 = mesh->indices[t * 3 + 2];

                int shared = 0;
                if (ti0 == i0 || ti0 == i1 || ti0 == i2) shared++;
                if (ti1 == i0 || ti1 == i1 || ti1 == i2) shared++;
                if (ti2 == i0 || ti2 == i1 || ti2 == i2) shared++;

                if (shared >= 2) {
                    /* Check normal similarity */
                    float dot = geomlib_vec3_dot(chart_normal, tri_normals[t]);
                    if (dot >= cos_threshold) {
                        assigned[t] = true;
                        queue[queue_size++] = t;
                    }
                }
            }
        }

        free(queue);
        atlas->chart_count++;

        if (atlas->chart_count >= 256) break;  /* Reasonable chart limit */
    }

    free(assigned);
    free(tri_normals);

    return GEOMLIB_MESH_SUCCESS;
}

/* Parameterize a single chart */
static int parameterize_chart(
    const geomlib_trimesh *mesh,
    geomlib_atlas_chart *chart,
    geomlib_param_method method
) {
    /* Create sub-mesh for this chart */
    geomlib_trimesh chart_mesh;

    /* Collect unique vertices */
    bool *vertex_used = calloc(mesh->vertex_count, sizeof(bool));
    for (size_t t = 0; t < chart->triangle_count; t++) {
        uint32_t tri_idx = chart->triangles[t];
        vertex_used[mesh->indices[tri_idx * 3 + 0]] = true;
        vertex_used[mesh->indices[tri_idx * 3 + 1]] = true;
        vertex_used[mesh->indices[tri_idx * 3 + 2]] = true;
    }

    /* Count and build vertex map */
    size_t vertex_count = 0;
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        if (vertex_used[v]) vertex_count++;
    }

    chart->vertex_count = vertex_count;
    chart->vertex_map = malloc(vertex_count * sizeof(uint32_t));
    uint32_t *vertex_remap = malloc(mesh->vertex_count * sizeof(uint32_t));

    size_t v_idx = 0;
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        if (vertex_used[v]) {
            chart->vertex_map[v_idx] = v;
            vertex_remap[v] = v_idx;
            v_idx++;
        }
    }

    /* Build chart mesh */
    geomlib_trimesh_init(&chart_mesh, vertex_count, chart->triangle_count);

    for (size_t i = 0; i < vertex_count; i++) {
        uint32_t orig_v = chart->vertex_map[i];
        geomlib_trimesh_add_vertex(&chart_mesh, mesh->vertices[orig_v], NULL, NULL);
    }

    for (size_t t = 0; t < chart->triangle_count; t++) {
        uint32_t tri_idx = chart->triangles[t];
        uint32_t i0 = vertex_remap[mesh->indices[tri_idx * 3 + 0]];
        uint32_t i1 = vertex_remap[mesh->indices[tri_idx * 3 + 1]];
        uint32_t i2 = vertex_remap[mesh->indices[tri_idx * 3 + 2]];
        geomlib_trimesh_add_triangle(&chart_mesh, i0, i1, i2);
    }

    /* Parameterize */
    geomlib_param_params uv_params = geomlib_param_params_default(method);
    geomlib_uv_coords uv_coords;
    uv_coords.uvs = NULL;
    uv_coords.vertex_count = 0;
    uv_coords.distortion = NULL;
    uv_coords.triangle_count = 0;

    int result = geomlib_parameterize(&chart_mesh, &uv_params, &uv_coords);

    if (result == GEOMLIB_MESH_SUCCESS && uv_coords.uvs) {
        /* Copy UVs */
        chart->uvs = malloc(vertex_count * sizeof(geomlib_vec2));
        memcpy(chart->uvs, uv_coords.uvs, vertex_count * sizeof(geomlib_vec2));

        /* Compute bounding box */
        chart->uv_min = (geomlib_vec2){FLT_MAX, FLT_MAX};
        chart->uv_max = (geomlib_vec2){-FLT_MAX, -FLT_MAX};

        for (size_t i = 0; i < vertex_count; i++) {
            geomlib_vec2 uv = chart->uvs[i];
            if (uv.x < chart->uv_min.x) chart->uv_min.x = uv.x;
            if (uv.y < chart->uv_min.y) chart->uv_min.y = uv.y;
            if (uv.x > chart->uv_max.x) chart->uv_max.x = uv.x;
            if (uv.y > chart->uv_max.y) chart->uv_max.y = uv.y;
        }
    }

    free(vertex_used);
    free(vertex_remap);
    geomlib_trimesh_free(&chart_mesh);
    if (uv_coords.uvs) free(uv_coords.uvs);
    if (uv_coords.distortion) free(uv_coords.distortion);

    return result;
}

/* Segment mesh into charts */
int geomlib_atlas_segment_mesh(
    const geomlib_trimesh *mesh,
    const geomlib_atlas_params *params,
    geomlib_atlas *atlas
) {
    if (!mesh || !params || !atlas) return GEOMLIB_MESH_ERROR_INVALID;

    switch (params->segmentation) {
        case GEOMLIB_ATLAS_SEG_NORMAL_CLUSTER:
        case GEOMLIB_ATLAS_SEG_AUTO:
            return segment_by_normals(mesh, params, atlas);

        default:
            /* Fallback: single chart */
            atlas->chart_count = 1;
            atlas->charts[0].triangles = malloc(mesh->triangle_count * sizeof(uint32_t));
            atlas->charts[0].triangle_count = mesh->triangle_count;
            for (size_t t = 0; t < mesh->triangle_count; t++) {
                atlas->charts[0].triangles[t] = t;
            }
            return GEOMLIB_MESH_SUCCESS;
    }
}

/* Parameterize all charts */
int geomlib_atlas_parameterize_charts(
    const geomlib_trimesh *mesh,
    geomlib_atlas *atlas,
    geomlib_param_method method
) {
    if (!mesh || !atlas) return GEOMLIB_MESH_ERROR_INVALID;

    for (size_t c = 0; c < atlas->chart_count; c++) {
        int result = parameterize_chart(mesh, &atlas->charts[c], method);
        if (result != GEOMLIB_MESH_SUCCESS) {
            return result;
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* Simple shelf packing algorithm */
static int pack_charts_shelf(
    geomlib_atlas *atlas,
    const geomlib_atlas_params *params
) {
    if (!atlas) return GEOMLIB_MESH_ERROR_INVALID;

    /* Sort charts by height (descending) for better packing */
    for (size_t i = 0; i < atlas->chart_count - 1; i++) {
        for (size_t j = i + 1; j < atlas->chart_count; j++) {
            float h_i = atlas->charts[i].uv_max.y - atlas->charts[i].uv_min.y;
            float h_j = atlas->charts[j].uv_max.y - atlas->charts[j].uv_min.y;
            if (h_j > h_i) {
                geomlib_atlas_chart temp = atlas->charts[i];
                atlas->charts[i] = atlas->charts[j];
                atlas->charts[j] = temp;
            }
        }
    }

    /* Shelf packing */
    float current_x = 0.0f;
    float current_y = 0.0f;
    float shelf_height = 0.0f;
    float gutter = params->gutter_width;

    for (size_t c = 0; c < atlas->chart_count; c++) {
        geomlib_atlas_chart *chart = &atlas->charts[c];

        /* Normalize UV to [0,1] */
        float width = chart->uv_max.x - chart->uv_min.x;
        float height = chart->uv_max.y - chart->uv_min.y;

        /* Check if we need a new shelf */
        if (current_x + width + gutter > 1.0f) {
            current_x = 0.0f;
            current_y += shelf_height + gutter;
            shelf_height = 0.0f;
        }

        /* Place chart */
        chart->atlas_offset = (geomlib_vec2){current_x, current_y};
        chart->atlas_scale = (geomlib_vec2){1.0f, 1.0f};

        /* Update position */
        current_x += width + gutter;
        if (height > shelf_height) {
            shelf_height = height;
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* Pack charts into atlas */
int geomlib_atlas_pack_charts(
    geomlib_atlas *atlas,
    const geomlib_atlas_params *params
) {
    if (!atlas || !params) return GEOMLIB_MESH_ERROR_INVALID;

    switch (params->packing) {
        case GEOMLIB_ATLAS_PACK_SHELF:
        default:
            return pack_charts_shelf(atlas, params);
    }
}

/* Apply atlas UVs to mesh */
int geomlib_atlas_apply_to_mesh(
    geomlib_trimesh *mesh,
    const geomlib_atlas *atlas
) {
    if (!mesh || !atlas) return GEOMLIB_MESH_ERROR_INVALID;

    /* Allocate UV array */
    if (!mesh->texcoords) {
        mesh->texcoords = malloc(mesh->vertex_count * sizeof(geomlib_vec2));
    }

    /* Apply each chart's UVs */
    for (size_t c = 0; c < atlas->chart_count; c++) {
        const geomlib_atlas_chart *chart = &atlas->charts[c];

        for (size_t i = 0; i < chart->vertex_count; i++) {
            uint32_t mesh_v = chart->vertex_map[i];
            geomlib_vec2 uv = chart->uvs[i];

            /* Transform UV to atlas space */
            uv.x = (uv.x - chart->uv_min.x) * chart->atlas_scale.x + chart->atlas_offset.x;
            uv.y = (uv.y - chart->uv_min.y) * chart->atlas_scale.y + chart->atlas_offset.y;

            mesh->texcoords[mesh_v] = uv;
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* Compute atlas statistics */
void geomlib_atlas_compute_stats(
    const geomlib_trimesh *mesh,
    geomlib_atlas *atlas
) {
    if (!mesh || !atlas) return;

    /* Compute coverage */
    float used_area = 0.0f;
    for (size_t c = 0; c < atlas->chart_count; c++) {
        const geomlib_atlas_chart *chart = &atlas->charts[c];
        float w = chart->uv_max.x - chart->uv_min.x;
        float h = chart->uv_max.y - chart->uv_min.y;
        used_area += w * h;
    }
    atlas->coverage = used_area * 100.0f;

    /* Compute stretch (simplified) */
    atlas->avg_stretch = 1.0f;
    atlas->max_stretch = 1.0f;
}

/* Get chart for triangle */
int geomlib_atlas_get_chart_for_triangle(
    const geomlib_atlas *atlas,
    uint32_t triangle_idx
) {
    if (!atlas) return -1;

    for (size_t c = 0; c < atlas->chart_count; c++) {
        for (size_t t = 0; t < atlas->charts[c].triangle_count; t++) {
            if (atlas->charts[c].triangles[t] == triangle_idx) {
                return (int)c;
            }
        }
    }

    return -1;
}

/* Generate complete atlas */
int geomlib_generate_atlas(
    geomlib_trimesh *mesh,
    const geomlib_atlas_params *params,
    geomlib_atlas *atlas
) {
    if (!mesh || !params || !atlas) return GEOMLIB_MESH_ERROR_INVALID;

    int result;

    /* Step 1: Segment mesh into charts */
    result = geomlib_atlas_segment_mesh(mesh, params, atlas);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Step 2: Parameterize each chart */
    result = geomlib_atlas_parameterize_charts(mesh, atlas, params->parameterization);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Step 3: Pack charts into atlas */
    result = geomlib_atlas_pack_charts(atlas, params);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Step 4: Apply UVs to mesh */
    result = geomlib_atlas_apply_to_mesh(mesh, atlas);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Step 5: Compute statistics */
    geomlib_atlas_compute_stats(mesh, atlas);

    return GEOMLIB_MESH_SUCCESS;
}
