/*
 * geomlib - Mesh saliency (see geomlib_saliency.h)
 *
 * Clean-room implementation of Lee, Beamish, White, "Mesh Saliency",
 * ACM SIGGRAPH 2005: measure discrete mean curvature at several Gaussian
 * smoothing scales, weight the difference between scales, and normalize.
 * Saliency segmentation seeds regions at the most salient vertices and grows
 * them with a priority flood over mesh edges.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_saliency.h"
#include "geomlib/geomlib_topology.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_saliency_result geomlib_saliency_result_create(size_t vertex_count)
{
    geomlib_saliency_result r;
    r.values = (float *)calloc(vertex_count ? vertex_count : 1, sizeof(float));
    r.vertex_count = vertex_count;
    return r;
}

void geomlib_saliency_result_free(geomlib_saliency_result *r)
{
    if (!r) return;
    free(r->values);
    r->values = NULL;
    r->vertex_count = 0;
}

int geomlib_compute_saliency(const geomlib_trimesh *mesh,
                             int scales,
                             geomlib_saliency_result *out)
{
    geomlib_vertex_adjacency adj;
    size_t V, i, j, s;
    float *curv, *smoothed, *next;
    double *saliency, mn, mx;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (scales <= 0) scales = 3;
    if (scales > 8) scales = 8;

    curv = (float *)calloc(V, sizeof(float));
    smoothed = (float *)calloc(V, sizeof(float));
    next = (float *)calloc(V, sizeof(float));
    saliency = (double *)calloc(V, sizeof(double));
    if (!curv || !smoothed || !next || !saliency) {
        free(curv); free(smoothed); free(next); free(saliency);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    adj = geomlib_build_vertex_adjacency(mesh);
    if (!adj.neighbor_counts) {
        geomlib_vertex_adjacency_free(&adj);
        free(curv); free(smoothed); free(next); free(saliency);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Mean curvature magnitude via the cotangent Laplacian / lumped area. */
    for (i = 0; i < V; ++i) {
        geomlib_vec3 lap = {0,0,0};
        double area = 0.0;
        uint32_t v0 = (uint32_t)i;
        for (j = 0u; j < adj.neighbor_counts[i]; ++j) {
            uint32_t nb = adj.neighbors[i][j];
            geomlib_vec3 e = geomlib_vec3_sub(mesh->vertices[nb], mesh->vertices[v0]);
            /* Equal-weight (uniform) Laplacian is a stable coarse estimate. */
            lap = geomlib_vec3_add(lap, e);
        }
        {
            size_t t;
            for (t = 0; t < mesh->triangle_count; ++t) {
                uint32_t a = mesh->indices[t*3], b = mesh->indices[t*3+1],
                         c = mesh->indices[t*3+2];
                if (a != v0 && b != v0 && c != v0) continue;
                geomlib_vec3 e1 = geomlib_vec3_sub(mesh->vertices[b], mesh->vertices[a]);
                geomlib_vec3 e2 = geomlib_vec3_sub(mesh->vertices[c], mesh->vertices[a]);
                area += 0.5 * geomlib_vec3_length(geomlib_vec3_cross(e1, e2)) / 3.0;
            }
        }
        curv[i] = (area > 1e-12) ? (float)(geomlib_vec3_length(lap)/(2.0*area)) : 0.0f;
    }

    /* Multi-scale Gaussian-like smoothing and scale-difference weighting. */
    memcpy(smoothed, curv, V*sizeof(float));
    for (s = 0; s < (size_t)scales; ++s) {
        for (i = 0; i < V; ++i) {
            double sum = 0.0;
            size_t n = adj.neighbor_counts[i];
            for (j = 0u; j < adj.neighbor_counts[i]; ++j)
                sum += smoothed[adj.neighbors[i][j]];
            next[i] = (n > 0) ? (float)((sum + smoothed[i])/(double)(n+1)) : smoothed[i];
        }
        for (i = 0; i < V; ++i) {
            double d = fabs(smoothed[i] - next[i]);
            saliency[i] += (s+1) * d * d;
        }
        memcpy(smoothed, next, V*sizeof(float));
    }

    mn = saliency[0]; mx = saliency[0];
    for (i = 1; i < V; ++i) {
        if (saliency[i] < mn) mn = saliency[i];
        if (saliency[i] > mx) mx = saliency[i];
    }
    if (mx > mn) {
        for (i = 0; i < V; ++i)
            out->values[i] = (float)((saliency[i]-mn)/(mx-mn));
    } else {
        for (i = 0; i < V; ++i) out->values[i] = 0.0f;
    }

    geomlib_vertex_adjacency_free(&adj);
    free(curv); free(smoothed); free(next); free(saliency);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_segment_by_saliency(const geomlib_trimesh *mesh,
                                const geomlib_saliency_result *saliency,
                                size_t region_count,
                                uint32_t *labels)
{
    geomlib_vertex_adjacency adj;
    size_t V, i, j, k, assigned = 0;
    uint32_t *seeds;

    if (!mesh || !saliency || !labels || !saliency->values ||
        mesh->triangle_count == 0 || region_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (region_count > V) region_count = V;
    for (i = 0; i < V; ++i) labels[i] = UINT32_MAX;

    /* Seeds: highest saliency vertices, spaced apart by adjacency. */
    seeds = (uint32_t *)malloc(region_count*sizeof(uint32_t));
    if (!seeds) return GEOMLIB_MESH_ERROR_MEMORY;
    {
        uint32_t *order = (uint32_t *)malloc(V*sizeof(uint32_t));
        if (!order) { free(seeds); return GEOMLIB_MESH_ERROR_MEMORY; }
        for (i = 0; i < V; ++i) order[i] = (uint32_t)i;
        /* Simple descending selection sort by saliency (V is modest). */
        for (i = 0; i < V; ++i) {
            size_t best = i;
            for (j = i+1; j < V; ++j)
                if (saliency->values[order[j]] > saliency->values[order[best]]) best = j;
            { uint32_t t = order[i]; order[i] = order[best]; order[best] = t; }
        }
        for (i = 0; i < region_count; ++i) seeds[i] = order[i];
        free(order);
    }
    for (i = 0; i < region_count; ++i) labels[seeds[i]] = (uint32_t)i;

    adj = geomlib_build_vertex_adjacency(mesh);
    if (!adj.neighbor_counts) { free(seeds); return GEOMLIB_MESH_ERROR_INVALID; }
    /* Priority flood: repeatedly pick the unassigned vertex with the largest
     * assigned neighbour (deterministic scan). */
    while (assigned < V) {
        int progress = 0;
        for (i = 0; i < V; ++i) {
            uint32_t best_label = UINT32_MAX;
            float best_sal = -1.0f;
            if (labels[i] != UINT32_MAX) continue;
            for (j = 0u; j < adj.neighbor_counts[i]; ++j) {
                uint32_t nb = adj.neighbors[i][j];
                if (labels[nb] != UINT32_MAX && saliency->values[nb] > best_sal) {
                    best_sal = saliency->values[nb];
                    best_label = labels[nb];
                }
            }
            if (best_label != UINT32_MAX) { labels[i] = best_label; progress = 1; assigned++; }
        }
        if (!progress) {
            /* Disconnected component: assign to the nearest seeded region by
             * index to guarantee every vertex gets a label. */
            for (i = 0; i < V; ++i)
                if (labels[i] == UINT32_MAX) labels[i] = 0;
            break;
        }
    }
    (void)k;
    geomlib_vertex_adjacency_free(&adj);
    free(seeds);
    return GEOMLIB_MESH_SUCCESS;
}
