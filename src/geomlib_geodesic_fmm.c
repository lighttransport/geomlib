/*
 * geomlib - Fast-marching-style surface geodesics (see header)
 *
 * Clean-room implementation: Dijkstra on the edge graph gives an upper bound,
 * then repeated triangle-interior relaxation approximates propagation across
 * triangle interiors. This is an approximation of the fast marching method on
 * triangulated domains (Sethian; Kimmel & Sethian 1998), using edge sampling
 * instead of an exact local solver.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_fmm.h"
#include "geomlib/geomlib_topology.h"

#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_geodesic_fmm_result geomlib_geodesic_fmm_result_create(size_t count)
{
    geomlib_geodesic_fmm_result r;
    r.distances = (float *)malloc((count?count:1)*sizeof(float));
    r.count = count;
    if (r.distances)
        for (size_t i = 0; i < count; ++i) r.distances[i] = INFINITY;
    return r;
}

void geomlib_geodesic_fmm_result_free(geomlib_geodesic_fmm_result *r)
{
    if (!r) return;
    free(r->distances);
    r->distances = NULL;
    r->count = 0;
}

int geomlib_geodesic_fast_marching(const geomlib_trimesh *mesh,
                                   uint32_t source,
                                   int relaxation_passes,
                                   geomlib_geodesic_fmm_result *out)
{
    geomlib_vertex_adjacency adj = {0};
    float *dist;
    bool *visited;
    size_t V, i, j, pass, s, k;
    const size_t samples = 16;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        source >= mesh->vertex_count || mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (relaxation_passes < 0) relaxation_passes = 0;
    if (relaxation_passes > 20) relaxation_passes = 20;

    dist = (float *)malloc(V*sizeof(float));
    visited = (bool *)calloc(V, sizeof(bool));
    if (!dist || !visited) { free(dist); free(visited); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < V; ++i) dist[i] = INFINITY;
    dist[source] = 0.0f;
    adj = geomlib_build_vertex_adjacency(mesh);
    if (!adj.neighbor_counts) {
        free(dist); free(visited);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Dijkstra on the edge graph. */
    for (i = 0; i < V; ++i) {
        size_t best = V;
        float bestd = INFINITY;
        for (j = 0; j < V; ++j)
            if (!visited[j] && dist[j] < bestd) { bestd = dist[j]; best = j; }
        if (best == V || !isfinite(bestd)) break;
        visited[best] = true;
        for (j = 0; j < adj.neighbor_counts[best]; ++j) {
            uint32_t nb = adj.neighbors[best][j];
            geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[nb], mesh->vertices[best]);
            float l = geomlib_vec3_length(d);
            if (dist[best] + l < dist[nb]) dist[nb] = dist[best] + l;
        }
    }

    /* Triangle-interior relaxation: a point on an edge carries the linearly
     * interpolated distance of its endpoints; a vertex receives the minimum
     * interpolated distance plus its in-triangle distance to that point. */
    for (pass = 0; pass < (size_t)relaxation_passes; ++pass) {
        for (s = 0; s < mesh->triangle_count; ++s) {
            uint32_t v[3] = { mesh->indices[s*3], mesh->indices[s*3+1],
                              mesh->indices[s*3+2] };
            geomlib_vec3 p[3];
            for (k = 0; k < 3; ++k) p[k] = mesh->vertices[v[k]];
            for (k = 0; k < 3; ++k) {
                uint32_t a = v[k], b = v[(k+1)%3], c = v[(k+2)%3];
                geomlib_vec3 pa = p[k], pb = p[(k+1)%3], pc = p[(k+2)%3];
                if (!isfinite(dist[a]) || !isfinite(dist[b])) continue;
                for (j = 1; j < samples; ++j) {
                    float t = (float)j/(float)samples;
                    geomlib_vec3 q = geomlib_vec3_add(geomlib_vec3_mul(pa,1.0f-t),
                                                      geomlib_vec3_mul(pb,t));
                    float dq = dist[a]*(1.0f-t) + dist[b]*t;
                    float dc = geomlib_vec3_length(geomlib_vec3_sub(pc, q));
                    if (dq + dc < dist[c]) dist[c] = dq + dc;
                }
            }
        }
    }

    geomlib_vertex_adjacency_free(&adj);
    {
        geomlib_geodesic_fmm_result result = geomlib_geodesic_fmm_result_create(V);
        if (!result.distances) { free(dist); free(visited); return GEOMLIB_MESH_ERROR_MEMORY; }
        memcpy(result.distances, dist, V*sizeof(float));
        free(out->distances);
        *out = result;
    }
    free(dist);
    free(visited);
    return GEOMLIB_MESH_SUCCESS;
}
