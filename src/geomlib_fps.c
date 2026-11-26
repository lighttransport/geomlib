/*
 * geomlib - Farthest-point sampling (see geomlib_fps.h)
 *
 * Clean-room greedy farthest-point traversal: at each step compute the
 * distance from every surface vertex to the nearest sample and select the
 * farthest one. Deterministic ties resolve to the lowest vertex index.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_fps.h"
#include "geomlib/geomlib_geodesic.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_fps_result geomlib_fps_result_create(void)
{
    geomlib_fps_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_fps_result_free(geomlib_fps_result *r)
{
    if (!r) return;
    free(r->indices);
    free(r->radii);
    memset(r, 0, sizeof(*r));
}

/* Geodesic distance to the nearest sample, with a per-source fallback. */
static void nearest_sample_distances(const geomlib_trimesh *mesh,
                                     const uint32_t *sources, size_t nsources,
                                     float *dist)
{
    geomlib_geodesic_distances g;
    size_t i;
    int rc;

    g = geomlib_geodesic_create(mesh->vertex_count);
    rc = geomlib_compute_geodesic_heat(mesh, sources, nsources, 0.0f, &g);
    if (rc == GEOMLIB_MESH_SUCCESS && g.distances) {
        memcpy(dist, g.distances, mesh->vertex_count*sizeof(float));
        geomlib_geodesic_free(&g);
        return;
    }
    geomlib_geodesic_free(&g);
    for (i = 0; i < mesh->vertex_count; ++i) dist[i] = FLT_MAX;
    for (i = 0; i < nsources; ++i) {
        g = geomlib_geodesic_create(mesh->vertex_count);
        if (geomlib_compute_geodesic_from_vertex(mesh, sources[i], &g) ==
            GEOMLIB_MESH_SUCCESS && g.distances) {
            size_t v;
            for (v = 0; v < mesh->vertex_count; ++v)
                if (g.distances[v] < dist[v]) dist[v] = g.distances[v];
        }
        geomlib_geodesic_free(&g);
    }
    for (i = 0; i < mesh->vertex_count; ++i)
        if (!isfinite(dist[i])) dist[i] = FLT_MAX;
}

int geomlib_farthest_point_sampling(const geomlib_trimesh *mesh, size_t k,
                                    uint32_t seed, geomlib_fps_result *out)
{
    geomlib_fps_result r = geomlib_fps_result_create();
    uint32_t *sources = NULL;
    float *dist = NULL;
    size_t i, best = 0;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->vertex_count || !k || seed >= mesh->vertex_count) return rc;
    if (k > mesh->vertex_count) k = mesh->vertex_count;
    r.indices = (uint32_t *)malloc(k*sizeof(uint32_t));
    r.radii = (float *)malloc(k*sizeof(float));
    sources = (uint32_t *)malloc(k*sizeof(uint32_t));
    dist = (float *)malloc(mesh->vertex_count*sizeof(float));
    if (!r.indices || !r.radii || !sources || !dist) goto cleanup;

    sources[0] = seed;
    r.indices[0] = seed;
    r.radii[0] = 0.0f;
    r.count = 1;
    for (i = 1; i < k; ++i) {
        float mx = -1.0f;
        size_t v;
        nearest_sample_distances(mesh, sources, i, dist);
        for (v = 0; v < mesh->vertex_count; ++v)
            if (dist[v] > mx) { mx = dist[v]; best = v; }
        if (!(mx > 1e-7f)) break;
        sources[i] = (uint32_t)best;
        r.indices[i] = (uint32_t)best;
        r.radii[i] = mx;
        r.count++;
    }
    /* The last recorded value was the radius used to select the final sample;
     * replace it with the coverage after that sample was added. */
    if (r.count > 1) {
        float mx = -1.0f;
        size_t v;
        nearest_sample_distances(mesh, sources, r.count, dist);
        for (v = 0; v < mesh->vertex_count; ++v)
            if (dist[v] > mx) mx = dist[v];
        if (mx >= 0.0f) r.radii[r.count-1] = mx;
    }
    *out = r;
    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    free(sources);
    free(dist);
    if (rc != GEOMLIB_MESH_SUCCESS) geomlib_fps_result_free(&r);
    return rc;
}
