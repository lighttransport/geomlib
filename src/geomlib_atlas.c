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
    geomlib_uv_coords uv_coords = geomlib_uv_coords_create(vertex_count,
                                                           chart->triangle_count);

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
    geomlib_uv_coords_free(&uv_coords);

    return result;
}

/* Generic region growing driven by a per-triangle acceptance predicate. */
typedef int (*chart_predicate)(const geomlib_trimesh *mesh, size_t seed,
                               size_t candidate, const void *context);

static int segment_by_predicate(const geomlib_trimesh *mesh,
                                geomlib_atlas *atlas,
                                chart_predicate accepts, const void *context)
{
    bool *assigned;
    size_t t, i;

    if (!mesh || !atlas || !accepts) return GEOMLIB_MESH_ERROR_INVALID;
    assigned = (bool *)calloc(mesh->triangle_count, sizeof(bool));
    if (!assigned) return GEOMLIB_MESH_ERROR_MEMORY;

    atlas->chart_count = 0;
    for (t = 0; t < mesh->triangle_count; ++t) {
        if (assigned[t]) continue;
        if (atlas->chart_count >= 256) break;

        geomlib_atlas_chart *chart = &atlas->charts[atlas->chart_count];
        chart->triangles = (uint32_t *)malloc(mesh->triangle_count * sizeof(uint32_t));
        chart->triangle_count = 0;
        if (!chart->triangles) { free(assigned); return GEOMLIB_MESH_ERROR_MEMORY; }

        uint32_t *queue = (uint32_t *)malloc(mesh->triangle_count * sizeof(uint32_t));
        if (!queue) { free(chart->triangles); free(assigned); return GEOMLIB_MESH_ERROR_MEMORY; }
        size_t qsize = 0;
        queue[qsize++] = (uint32_t)t;
        assigned[t] = true;

        while (qsize > 0) {
            uint32_t current = queue[--qsize];
            chart->triangles[chart->triangle_count++] = current;
            for (i = 0; i < mesh->triangle_count; ++i) {
                if (assigned[i]) continue;
                if (accepts(mesh, current, i, context)) {
                    assigned[i] = true;
                    queue[qsize++] = (uint32_t)i;
                }
            }
        }
        free(queue);
        atlas->chart_count++;
    }

    free(assigned);
    return GEOMLIB_MESH_SUCCESS;
}

/* Planarity criterion: similar normal and vertices close to the chart plane. */
typedef struct {
    geomlib_vec3 *normals;
    float cos_threshold;
    float distance_epsilon;
} planarity_context;

static int accepts_planar(const geomlib_trimesh *mesh, size_t seed,
                          size_t candidate, const void *ctx)
{
    const planarity_context *c = (const planarity_context *)ctx;
    uint32_t a = mesh->indices[seed*3], b = mesh->indices[seed*3+1], cc = mesh->indices[seed*3+2];
    geomlib_vec3 point = mesh->vertices[a];
    geomlib_vec3 normal = c->normals[seed];
    (void)b; (void)cc;
    if (geomlib_vec3_dot(normal, c->normals[candidate]) < c->cos_threshold) return 0;
    for (int k = 0; k < 3; ++k) {
        uint32_t v = mesh->indices[candidate*3+k];
        geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[v], point);
        if (fabsf(geomlib_vec3_dot(d, normal)) > c->distance_epsilon) return 0;
    }
    return 1;
}

static int segment_by_planarity(const geomlib_trimesh *mesh,
                                const geomlib_atlas_params *params,
                                geomlib_atlas *atlas)
{
    planarity_context ctx;
    size_t t;
    ctx.normals = (geomlib_vec3 *)malloc(mesh->triangle_count * sizeof(geomlib_vec3));
    if (!ctx.normals) return GEOMLIB_MESH_ERROR_MEMORY;
    for (t = 0; t < mesh->triangle_count; ++t) {
        uint32_t i0 = mesh->indices[t*3], i1 = mesh->indices[t*3+1], i2 = mesh->indices[t*3+2];
        ctx.normals[t] = compute_triangle_normal(mesh->vertices[i0], mesh->vertices[i1],
                                                 mesh->vertices[i2]);
    }
    ctx.cos_threshold = cosf(params->normal_threshold * (float)M_PI / 180.0f);
    {
        geomlib_aabb box = geomlib_trimesh_compute_aabb(mesh);
        geomlib_vec3 size = geomlib_aabb_size(box);
        float diag = geomlib_vec3_length(size);
        ctx.distance_epsilon = 0.01f * diag;
    }
    int rc = segment_by_predicate(mesh, atlas, accepts_planar, &ctx);
    free(ctx.normals);
    return rc;
}

/* Distortion criterion: triangles with comparable shape anisotropy. */
static float triangle_anisotropy(const geomlib_trimesh *mesh, size_t t)
{
    uint32_t i0 = mesh->indices[t*3], i1 = mesh->indices[t*3+1], i2 = mesh->indices[t*3+2];
    float a = geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[i1], mesh->vertices[i0]));
    float b = geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[i2], mesh->vertices[i1]));
    float c = geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[i0], mesh->vertices[i2]));
    float lo = a < b ? a : b; if (c < lo) lo = c;
    float hi = a > b ? a : b; if (c > hi) hi = c;
    return (lo > 1e-8f) ? hi / lo : 1.0f;
}

static int accepts_similar_distortion(const geomlib_trimesh *mesh, size_t seed,
                                      size_t candidate, const void *ctx)
{
    const geomlib_atlas_params *p = (const geomlib_atlas_params *)ctx;
    float a = triangle_anisotropy(mesh, seed);
    float b = triangle_anisotropy(mesh, candidate);
    float limit = 1.0f + 0.5f * (p->distortion_threshold - 1.0f);
    return (a <= b * limit) && (b <= a * limit);
}

static int segment_by_distortion(const geomlib_trimesh *mesh,
                                 const geomlib_atlas_params *params,
                                 geomlib_atlas *atlas)
{
    return segment_by_predicate(mesh, atlas, accepts_similar_distortion, params);
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
        case GEOMLIB_ATLAS_SEG_PLANARITY:
            return segment_by_planarity(mesh, params, atlas);
        case GEOMLIB_ATLAS_SEG_DISTORTION:
            return segment_by_distortion(mesh, params, atlas);

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

/* Rotate chart UVs to minimize bounding-box area (chart-local copy only). */
static void optimize_chart_rotation(geomlib_atlas_chart *chart)
{
    float best_area = FLT_MAX, best_angle = 0.0f;
    geomlib_vec2 best_min = {0,0}, best_max = {1,1};
    int step;

    if (!chart->uvs || chart->vertex_count == 0) return;
    for (step = 0; step < 12; ++step) {
        float angle = step * (float)M_PI / 12.0f;
        float cs = cosf(angle), sn = sinf(angle);
        float min_u = FLT_MAX, min_v = FLT_MAX, max_u = -FLT_MAX, max_v = -FLT_MAX;
        for (size_t i = 0; i < chart->vertex_count; ++i) {
            float u = chart->uvs[i].x - chart->uv_min.x;
            float v = chart->uvs[i].y - chart->uv_min.y;
            float ru = cs*u - sn*v;
            float rv = sn*u + cs*v;
            if (ru < min_u) min_u = ru;
            if (rv < min_v) min_v = rv;
            if (ru > max_u) max_u = ru;
            if (rv > max_v) max_v = rv;
        }
        float area = (max_u - min_u) * (max_v - min_v);
        if (area < best_area) {
            best_area = area;
            best_angle = angle;
            best_min = (geomlib_vec2){min_u, min_v};
            best_max = (geomlib_vec2){max_u, max_v};
        }
    }
    {
        float cs = cosf(best_angle), sn = sinf(best_angle);
        for (size_t i = 0; i < chart->vertex_count; ++i) {
            float u = chart->uvs[i].x - chart->uv_min.x;
            float v = chart->uvs[i].y - chart->uv_min.y;
            chart->uvs[i].x = cs*u - sn*v - best_min.x;
            chart->uvs[i].y = sn*u + cs*v - best_min.y;
        }
    }
    chart->uv_min = (geomlib_vec2){0, 0};
    chart->uv_max = (geomlib_vec2){best_max.x - best_min.x, best_max.y - best_min.y};
}

typedef struct { float x, y, w, h; } free_rect;

#define GEOMLIB_ATLAS_MAX_FREE 4096

static void split_free_rect(free_rect *free_rects, size_t *count,
                            const free_rect *used, float gutter)
{
    free_rect result[GEOMLIB_ATLAS_MAX_FREE];
    size_t n = 0, i;
    for (i = 0; i < *count && n < GEOMLIB_ATLAS_MAX_FREE; ++i) {
        free_rect f = free_rects[i];
        float fx2 = f.x + f.w, fy2 = f.y + f.h;
        float ux2 = used->x + used->w, uy2 = used->y + used->h;
        free_rect r;
        if (used->x >= fx2 || ux2 <= f.x || used->y >= fy2 || uy2 <= f.y) {
            result[n++] = f;
            continue;
        }
        if (used->x > f.x) {
            r = (free_rect){f.x, f.y, used->x - f.x - gutter, f.h};
            if (r.w > 0.0f && n < GEOMLIB_ATLAS_MAX_FREE) result[n++] = r;
        }
        if (ux2 < fx2) {
            r = (free_rect){ux2 + gutter, f.y, fx2 - ux2 - gutter, f.h};
            if (r.w > 0.0f && n < GEOMLIB_ATLAS_MAX_FREE) result[n++] = r;
        }
        if (used->y > f.y) {
            r = (free_rect){f.x, f.y, f.w, used->y - f.y - gutter};
            if (r.h > 0.0f && n < GEOMLIB_ATLAS_MAX_FREE) result[n++] = r;
        }
        if (uy2 < fy2) {
            r = (free_rect){f.x, uy2 + gutter, f.w, fy2 - uy2 - gutter};
            if (r.h > 0.0f && n < GEOMLIB_ATLAS_MAX_FREE) result[n++] = r;
        }
    }
    memcpy(free_rects, result, n * sizeof(free_rect));
    *count = n;
}

/* MaxRects best-short-side-fit packing. */
static int pack_charts_maxrects(geomlib_atlas *atlas,
                                const geomlib_atlas_params *params)
{
    free_rect *free_rects = (free_rect *)malloc(GEOMLIB_ATLAS_MAX_FREE * sizeof(free_rect));
    size_t count = 0;
    float gutter = params->gutter_width;

    if (!free_rects) return GEOMLIB_MESH_ERROR_MEMORY;
    free_rects[count++] = (free_rect){0, 0, 1.0f, 1.0f};

    for (size_t c = 0; c < atlas->chart_count; ++c) {
        geomlib_atlas_chart *chart = &atlas->charts[c];
        float w = chart->uv_max.x - chart->uv_min.x + gutter;
        float h = chart->uv_max.y - chart->uv_min.y + gutter;
        size_t best = count, i;
        float best_score = FLT_MAX;
        for (i = 0; i < count; ++i) {
            if (free_rects[i].w < w || free_rects[i].h < h) continue;
            float short_side = free_rects[i].w - w < free_rects[i].h - h
                ? free_rects[i].w - w : free_rects[i].h - h;
            if (short_side < best_score) { best_score = short_side; best = i; }
        }
        if (best == count) {
            /* No rectangle is large enough: shrink into the largest free one. */
            size_t largest = 0;
            for (i = 1; i < count; ++i) {
                if (free_rects[i].w * free_rects[i].h >
                    free_rects[largest].w * free_rects[largest].h) largest = i;
            }
            if (count == 0) { chart->atlas_offset = (geomlib_vec2){0,0}; continue; }
            float scale = free_rects[largest].w / w;
            float sy = free_rects[largest].h / h;
            if (sy < scale) scale = sy;
            if (scale > 1.0f) scale = 1.0f;
            chart->atlas_scale = (geomlib_vec2){scale, scale};
            chart->atlas_offset = (geomlib_vec2){free_rects[largest].x,
                                                  free_rects[largest].y};
            free_rect used = (free_rect){free_rects[largest].x, free_rects[largest].y,
                                         w * scale, h * scale};
            split_free_rect(free_rects, &count, &used, gutter);
            continue;
        }
        {
            free_rect used = (free_rect){free_rects[best].x, free_rects[best].y, w, h};
            chart->atlas_offset = (geomlib_vec2){used.x, used.y};
            chart->atlas_scale = (geomlib_vec2){1, 1};
            split_free_rect(free_rects, &count, &used, gutter);
        }
    }
    free(free_rects);
    return GEOMLIB_MESH_SUCCESS;
}

/* Guillotine packing: place into the first fitting free rectangle, then split. */
static int pack_charts_guillotine(geomlib_atlas *atlas,
                                  const geomlib_atlas_params *params)
{
    free_rect *free_rects = (free_rect *)malloc(GEOMLIB_ATLAS_MAX_FREE * sizeof(free_rect));
    size_t count = 0;
    float gutter = params->gutter_width;
    if (!free_rects) return GEOMLIB_MESH_ERROR_MEMORY;
    free_rects[count++] = (free_rect){0, 0, 1.0f, 1.0f};

    for (size_t c = 0; c < atlas->chart_count; ++c) {
        geomlib_atlas_chart *chart = &atlas->charts[c];
        float w = chart->uv_max.x - chart->uv_min.x + gutter;
        float h = chart->uv_max.y - chart->uv_min.y + gutter;
        size_t chosen = count, i;
        for (i = 0; i < count; ++i) {
            if (free_rects[i].w >= w && free_rects[i].h >= h) { chosen = i; break; }
        }
        if (chosen == count) {
            size_t largest = 0;
            for (i = 1; i < count; ++i) {
                if (free_rects[i].w * free_rects[i].h >
                    free_rects[largest].w * free_rects[largest].h) largest = i;
            }
            if (count == 0) { chart->atlas_offset = (geomlib_vec2){0,0}; continue; }
            float scale = free_rects[largest].w / w;
            float sy = free_rects[largest].h / h;
            if (sy < scale) scale = sy;
            if (scale > 1.0f) scale = 1.0f;
            chart->atlas_scale = (geomlib_vec2){scale, scale};
            chart->atlas_offset = (geomlib_vec2){free_rects[largest].x,
                                                  free_rects[largest].y};
            continue;
        }
        {
            free_rect f = free_rects[chosen];
            chart->atlas_offset = (geomlib_vec2){f.x, f.y};
            chart->atlas_scale = (geomlib_vec2){1, 1};
            /* Remove the chosen rectangle and append the two guillotine cuts. */
            for (i = chosen; i + 1 < count; ++i) free_rects[i] = free_rects[i+1];
            count--;
            if (f.w - w > gutter && count < GEOMLIB_ATLAS_MAX_FREE)
                free_rects[count++] = (free_rect){f.x + w + gutter, f.y, f.w - w - gutter, f.h};
            if (f.h - h > gutter && count < GEOMLIB_ATLAS_MAX_FREE)
                free_rects[count++] = (free_rect){f.x, f.y + h + gutter, f.w, f.h - h - gutter};
        }
    }
    free(free_rects);
    return GEOMLIB_MESH_SUCCESS;
}

/* Pack charts into atlas */
int geomlib_atlas_pack_charts(
    geomlib_atlas *atlas,
    const geomlib_atlas_params *params
) {
    if (!atlas || !params) return GEOMLIB_MESH_ERROR_INVALID;

    if (params->optimize_rotation) {
        for (size_t c = 0; c < atlas->chart_count; ++c)
            optimize_chart_rotation(&atlas->charts[c]);
    }
    switch (params->packing) {
        case GEOMLIB_ATLAS_PACK_MAXRECTS:
            return pack_charts_maxrects(atlas, params);
        case GEOMLIB_ATLAS_PACK_GUILLOTINE:
            return pack_charts_guillotine(atlas, params);
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

    /* Coverage = packed chart area / atlas bounding area, as a percentage. */
    float used_area = 0.0f, extent_x = 0.0f, extent_y = 0.0f;
    for (size_t c = 0; c < atlas->chart_count; c++) {
        const geomlib_atlas_chart *chart = &atlas->charts[c];
        float w = (chart->uv_max.x - chart->uv_min.x) * chart->atlas_scale.x;
        float h = (chart->uv_max.y - chart->uv_min.y) * chart->atlas_scale.y;
        used_area += w * h;
        float x2 = chart->atlas_offset.x + w;
        float y2 = chart->atlas_offset.y + h;
        if (x2 > extent_x) extent_x = x2;
        if (y2 > extent_y) extent_y = y2;
    }
    if (extent_x > 0.0f && extent_y > 0.0f)
        atlas->coverage = 100.0f * used_area / (extent_x * extent_y);
    else
        atlas->coverage = 0.0f;

    /* Compute stretch: normalized UV-to-3D edge-length ratios. */
    double total_3d_area = 0.0, total_uv_area = 0.0;
    for (size_t c = 0; c < atlas->chart_count; c++) {
        const geomlib_atlas_chart *chart = &atlas->charts[c];
        for (size_t t = 0; t < chart->triangle_count; t++) {
            size_t triangle = chart->triangles[t];
            if (triangle >= mesh->triangle_count) continue;
            uint32_t i0 = mesh->indices[triangle * 3], i1 = mesh->indices[triangle * 3 + 1],
                     i2 = mesh->indices[triangle * 3 + 2];
            geomlib_vec3 p0 = mesh->vertices[i0], p1 = mesh->vertices[i1],
                          p2 = mesh->vertices[i2];
            total_3d_area += 0.5 * geomlib_vec3_length(geomlib_vec3_cross(
                geomlib_vec3_sub(p1, p0), geomlib_vec3_sub(p2, p0)));
            if (!mesh->texcoords) continue;
            geomlib_vec2 uv0 = mesh->texcoords[i0], uv1 = mesh->texcoords[i1],
                          uv2 = mesh->texcoords[i2];
            total_uv_area += 0.5 * fabs((uv1.x - uv0.x) * (uv2.y - uv0.y) -
                                        (uv1.y - uv0.y) * (uv2.x - uv0.x));
        }
    }
    double scale = 1.0;
    if (total_3d_area > 0.0 && total_uv_area > 0.0)
        scale = sqrt(total_uv_area / total_3d_area);
    double stretch_sum = 0.0;
    float max_stretch = 0.0f;
    size_t stretch_samples = 0;
    if (mesh->texcoords) {
        for (size_t c = 0; c < atlas->chart_count; c++) {
            const geomlib_atlas_chart *chart = &atlas->charts[c];
            for (size_t t = 0; t < chart->triangle_count; t++) {
                size_t triangle = chart->triangles[t];
                if (triangle >= mesh->triangle_count) continue;
                uint32_t index[3] = {mesh->indices[triangle * 3],
                                     mesh->indices[triangle * 3 + 1],
                                     mesh->indices[triangle * 3 + 2]};
                for (int j = 0; j < 3; j++) {
                    uint32_t a = index[j], b = index[(j + 1) % 3];
                    double length3d = geomlib_vec3_length(
                        geomlib_vec3_sub(mesh->vertices[b], mesh->vertices[a]));
                    double length_uv = hypot((double)mesh->texcoords[b].x - mesh->texcoords[a].x,
                                             (double)mesh->texcoords[b].y - mesh->texcoords[a].y);
                    if (!(length3d > 1e-12) || !(scale > 1e-12)) continue;
                    double ratio = length_uv / (length3d * scale);
                    stretch_sum += ratio;
                    if (ratio > max_stretch) max_stretch = (float)ratio;
                    stretch_samples++;
                }
            }
        }
    }
    if (stretch_samples) {
        atlas->avg_stretch = (float)(stretch_sum / (double)stretch_samples);
        atlas->max_stretch = max_stretch;
    } else {
        atlas->avg_stretch = 1.0f;
        atlas->max_stretch = 1.0f;
    }
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
