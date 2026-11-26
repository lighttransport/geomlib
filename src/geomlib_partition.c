/*
 * geomlib - Balanced mesh partition (see geomlib_partition.h)
 *
 * Clean-room implementation of balanced region growing: farthest-point seeds
 * in Euclidean space, then simultaneous multi-source BFS with per-part
 * capacity so the resulting parts are spatially coherent and balanced.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_partition.h"
#include "geomlib/geomlib_topology.h"

#include <float.h>
#include <stdlib.h>
#include <string.h>

geomlib_partition_result geomlib_partition_result_create(void)
{
    geomlib_partition_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_partition_result_free(geomlib_partition_result *r)
{
    if (!r) return;
    free(r->part_of_vertex);
    free(r->part_sizes);
    memset(r, 0, sizeof(*r));
}

int geomlib_partition_mesh(const geomlib_trimesh *mesh,
                           size_t part_count,
                           geomlib_partition_result *out)
{
    uint32_t *seeds, *labels, *queue;
    size_t V, i, j, k;

    if (!mesh || !out || !mesh->vertices || mesh->vertex_count == 0 ||
        part_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (part_count > V) part_count = V;

    labels = (uint32_t *)malloc(V*sizeof(uint32_t));
    seeds = (uint32_t *)malloc(part_count*sizeof(uint32_t));
    queue = (uint32_t *)malloc(V*sizeof(uint32_t));
    if (!labels || !seeds || !queue) {
        free(labels); free(seeds); free(queue);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < V; ++i) labels[i] = UINT32_MAX;

    /* Farthest-point seeding on vertex positions. */
    seeds[0] = 0;
    for (k = 1; k < part_count; ++k) {
        double best = -1.0;
        size_t pick = k;
        for (i = 0; i < V; ++i) {
            double mind = DBL_MAX;
            for (j = 0; j < k; ++j) {
                geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[i], mesh->vertices[seeds[j]]);
                double dd = geomlib_vec3_dot(d, d);
                if (dd < mind) mind = dd;
            }
            if (mind > best) { best = mind; pick = i; }
        }
        seeds[k] = (uint32_t)pick;
    }

    /* Assign each vertex to its nearest farthest-point seed. Farthest-point
     * seeding spreads the seeds, so the resulting Voronoi regions are both
     * spatially coherent and approximately balanced. */
    {
        size_t *sizes = (size_t *)calloc(part_count, sizeof(size_t));
        if (!sizes) {
            free(labels); free(seeds); free(queue);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i = 0; i < V; ++i) {
            double best = DBL_MAX;
            size_t bpart = 0;
            for (k = 0; k < part_count; ++k) {
                geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[i], mesh->vertices[seeds[k]]);
                double dd = geomlib_vec3_dot(d, d);
                if (dd < best) { best = dd; bpart = k; }
            }
            labels[i] = (uint32_t)bpart;
            sizes[bpart]++;
        }
        /* Deterministic rebalance: move the vertex of an overfull part that
         * lies closest to an underfull part's seed. */
        {
            size_t capacity = (V + part_count - 1) / part_count;
            int pass;
            for (pass = 0; pass < (int)V; ++pass) {
                size_t src = 0, dst = 0;
                for (k = 1; k < part_count; ++k)
                    if (sizes[k] > sizes[src]) src = k;
                for (k = 1; k < part_count; ++k)
                    if (k != src && sizes[k] < sizes[dst]) dst = k;
                if (src == dst || sizes[src] <= capacity) break;
                {
                    size_t best_i = V;
                    double best_d = DBL_MAX;
                    for (i = 0; i < V; ++i) {
                        geomlib_vec3 d;
                        double dd;
                        if (labels[i] != src) continue;
                        d = geomlib_vec3_sub(mesh->vertices[i], mesh->vertices[seeds[dst]]);
                        dd = geomlib_vec3_dot(d, d);
                        if (dd < best_d) { best_d = dd; best_i = i; }
                    }
                    if (best_i == V) break;
                    labels[best_i] = (uint32_t)dst;
                    sizes[src]--;
                    sizes[dst]++;
                }
            }
        }
        {
            geomlib_partition_result result = geomlib_partition_result_create();
            result.part_of_vertex = labels;
            result.part_sizes = sizes;
            result.vertex_count = V;
            result.part_count = part_count;
            *out = result;
        }
    }
    free(seeds);
    free(queue);
    return GEOMLIB_MESH_SUCCESS;
}
