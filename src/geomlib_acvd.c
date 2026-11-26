/*
 * geomlib - ACVD-style remeshing (see geomlib_acvd.h)
 *
 * Clean-room implementation of the clustering/CVT idea from:
 *   Valette, Chassery, Prost, "Generic remeshing of 3D triangular meshes with
 *   metric-dependent discrete Voronoi diagrams", IEEE TVCG 2008.
 *
 * Steps: dense area-weighted sampling, Lloyd relaxation of sites toward the
 * centroids of their nearest-neighbor clusters (projected onto the input
 * surface via BVH), then local tangent-plane angular triangulation of the
 * sites. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_acvd.h"
#include "geomlib/geomlib_bvh.h"
#include "geomlib/geomlib_random.h"
#include "geomlib/geomlib_repair.h"
#include "geomlib/geomlib_topology.h"
#include "geomlib/geomlib_surface.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_acvd_params geomlib_acvd_params_default(uint32_t target_vertices)
{
    geomlib_acvd_params p;
    p.target_vertices = target_vertices;
    p.lloyd_iterations = 5;
    p.seed = 0x12345678u;
    return p;
}

int geomlib_remesh_acvd(const geomlib_trimesh *input,
                        const geomlib_acvd_params *params,
                        geomlib_trimesh *output)
{
    geomlib_bvh bvh;
    geomlib_pcg32 rng;
    double *cdf;
    double total_area = 0.0;
    size_t t, i, j, k;
    uint32_t nsites, nsamples;
    geomlib_vec3 *sites = NULL, *normals = NULL, *samples = NULL;
    size_t *assignment = NULL;
    double *cx = NULL, *cy = NULL, *cz = NULL;
    int *counts = NULL;
    double avg_area, spacing;
    int rc = GEOMLIB_MESH_ERROR_INVALID;
    int iters;

    if (!input || !output || !input->vertices || !input->indices ||
        input->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    nsites = params ? params->target_vertices : 0;
    if (nsites < 10) nsites = 10;
    if (nsites > 4000) nsites = 4000;
    iters = params ? params->lloyd_iterations : 5;
    if (iters < 0) iters = 0;
    if (iters > 20) iters = 20;

    cdf = (double *)malloc(input->triangle_count * sizeof(double));
    if (!cdf) return GEOMLIB_MESH_ERROR_MEMORY;
    for (t = 0; t < input->triangle_count; ++t) {
        uint32_t a = input->indices[t*3], b = input->indices[t*3+1], c = input->indices[t*3+2];
        geomlib_vec3 e1 = geomlib_vec3_sub(input->vertices[b], input->vertices[a]);
        geomlib_vec3 e2 = geomlib_vec3_sub(input->vertices[c], input->vertices[a]);
        total_area += 0.5 * geomlib_vec3_length(geomlib_vec3_cross(e1, e2));
        cdf[t] = total_area;
    }
    if (!(total_area > 0.0)) { free(cdf); return GEOMLIB_MESH_ERROR_INVALID; }
    avg_area = total_area / nsites;
    spacing = sqrt(4.0 * avg_area / sqrt(3.0));

    geomlib_pcg32_init(&rng, params ? params->seed : 0x12345678u, 0x9e3779b9u);
    geomlib_bvh_build(&bvh, input, GEOMLIB_BVH_QUALITY_MEDIUM);
    nsamples = nsites * 8;
    sites = (geomlib_vec3 *)malloc(nsites * sizeof(geomlib_vec3));
    normals = (geomlib_vec3 *)malloc(nsites * sizeof(geomlib_vec3));
    samples = (geomlib_vec3 *)malloc(nsamples * sizeof(geomlib_vec3));
    assignment = (size_t *)malloc(nsamples * sizeof(size_t));
    cx = (double *)calloc(nsites, sizeof(double));
    cy = (double *)calloc(nsites, sizeof(double));
    cz = (double *)calloc(nsites, sizeof(double));
    counts = (int *)calloc(nsites, sizeof(int));
    if (!sites || !normals || !samples || !assignment || !cx || !cy || !cz || !counts) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }

#define SAMPLE_POINT(out) do {                                                   \
        double pick = geomlib_pcg32_nextd(&rng) * total_area;                    \
        size_t lo = 0, hi = input->triangle_count;                               \
        while (lo < hi) {                                                        \
            size_t mid = (lo + hi) / 2;                                          \
            if (cdf[mid] < pick) lo = mid + 1; else hi = mid;                    \
        }                                                                        \
        if (lo >= input->triangle_count) lo = input->triangle_count - 1;         \
        {                                                                        \
            double u = geomlib_pcg32_nextd(&rng), v = geomlib_pcg32_nextd(&rng); \
            if (u + v > 1.0) { u = 1.0 - u; v = 1.0 - v; }                       \
            uint32_t a = input->indices[lo*3], b = input->indices[lo*3+1];       \
            uint32_t c = input->indices[lo*3+2];                                 \
            (out).x = input->vertices[a].x + (float)u*(input->vertices[b].x-input->vertices[a].x) \
                                           + (float)v*(input->vertices[c].x-input->vertices[a].x); \
            (out).y = input->vertices[a].y + (float)u*(input->vertices[b].y-input->vertices[a].y) \
                                           + (float)v*(input->vertices[c].y-input->vertices[a].y); \
            (out).z = input->vertices[a].z + (float)u*(input->vertices[b].z-input->vertices[a].z) \
                                           + (float)v*(input->vertices[c].z-input->vertices[a].z); \
        }                                                                        \
    } while (0)

    for (i = 0; i < nsites; ++i) {
        SAMPLE_POINT(sites[i]);
        {
            geomlib_vec3 closest;
            float dist;
            if (geomlib_bvh_closest_point(&bvh, sites[i], &closest, &dist, NULL))
                sites[i] = closest;
            normals[i] = geomlib_trimesh_triangle_normal(
                input, 0); /* replaced below with closest-face normal */
            /* Determine the closest face normal by finding the nearest vertex
               triangle center; adequate for tangent-plane construction. */
            {
                float best = 1e30f;
                size_t bt = 0;
                for (t = 0; t < input->triangle_count; ++t) {
                    uint32_t a = input->indices[t*3], b = input->indices[t*3+1], c = input->indices[t*3+2];
                    geomlib_vec3 fc = geomlib_vec3_mul(
                        geomlib_vec3_add(geomlib_vec3_add(input->vertices[a], input->vertices[b]),
                                         input->vertices[c]), 1.0f/3.0f);
                    geomlib_vec3 d = geomlib_vec3_sub(sites[i], fc);
                    float dd = geomlib_vec3_dot(d, d);
                    if (dd < best) { best = dd; bt = t; }
                }
                normals[i] = geomlib_trimesh_triangle_normal(input, bt);
            }
        }
    }
    for (i = 0; i < nsamples; ++i) SAMPLE_POINT(samples[i]);
#undef SAMPLE_POINT

    /* Lloyd relaxation: assign dense samples to nearest sites and move each
     * site to the centroid of its cluster, then project back to the surface. */
    for (i = 0; i < (size_t)iters; ++i) {
        memset(cx, 0, nsites*sizeof(double));
        memset(cy, 0, nsites*sizeof(double));
        memset(cz, 0, nsites*sizeof(double));
        memset(counts, 0, nsites*sizeof(int));
        for (j = 0; j < nsamples; ++j) {
            size_t best = 0;
            float best_d = 1e30f;
            for (k = 0; k < nsites; ++k) {
                geomlib_vec3 d = geomlib_vec3_sub(samples[j], sites[k]);
                float dd = geomlib_vec3_dot(d, d);
                if (dd < best_d) { best_d = dd; best = k; }
            }
            assignment[j] = best;
            cx[best] += samples[j].x;
            cy[best] += samples[j].y;
            cz[best] += samples[j].z;
            counts[best]++;
        }
        for (k = 0; k < nsites; ++k) {
            if (counts[k] == 0) continue;
            geomlib_vec3 centroid = {
                (float)(cx[k]/counts[k]), (float)(cy[k]/counts[k]), (float)(cz[k]/counts[k])
            };
            geomlib_vec3 closest;
            float dist;
            if (geomlib_bvh_closest_point(&bvh, centroid, &closest, &dist, NULL))
                sites[k] = closest;
        }
    }

    /* Triangulate the relaxed sites with the existing greedy projection
     * reconstruction (tangent-plane fans with normal/edge checks), which
     * produces a single consistent surface rather than overlapping fans. */
    {
        geomlib_pointcloud pc;
        geomlib_surface_params sp;
        int prc = geomlib_pointcloud_init(&pc, nsites);
        if (prc != GEOMLIB_MESH_SUCCESS) {
            rc = GEOMLIB_MESH_ERROR_MEMORY;
            goto cleanup;
        }
        for (k = 0; k < nsites; ++k)
            geomlib_pointcloud_add_point(&pc, sites[k], &normals[k], 1.0f);
        sp = geomlib_surface_params_default(GEOMLIB_SURFACE_GREEDY_PROJECTION);
        sp.search_radius = (float)(3.0 * spacing);
        sp.max_neighbors = 24;
        sp.max_edge_length = (float)(2.0 * spacing);
        sp.max_surface_angle = 60.0f;
        rc = geomlib_reconstruct_surface(&pc, &sp, output);
        geomlib_pointcloud_free(&pc);
    }

cleanup:
    free(cdf);
    free(sites); free(normals); free(samples); free(assignment);
    free(cx); free(cy); free(cz); free(counts);
    geomlib_bvh_free(&bvh);
    return rc;
}
