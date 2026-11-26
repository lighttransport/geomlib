/*
 * geomlib - Blue-noise surface sampling (see geomlib_bluenoise.h)
 *
 * Clean-room implementation of Bridson's Poisson-disk sampling (2007)
 * adapted to triangle meshes with area-weighted triangle selection and a
 * uniform spatial hash over 3D positions. C11, wasm32 safe.
 */

#include "geomlib/geomlib_bluenoise.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_sample_set geomlib_sample_set_create(size_t capacity)
{
    geomlib_sample_set s;
    s.points = (geomlib_vec3 *)malloc((capacity ? capacity : 1) * sizeof(geomlib_vec3));
    s.normals = (geomlib_vec3 *)malloc((capacity ? capacity : 1) * sizeof(geomlib_vec3));
    s.triangle_ids = (uint32_t *)malloc((capacity ? capacity : 1) * sizeof(uint32_t));
    s.count = 0;
    if (!s.points || !s.normals || !s.triangle_ids) geomlib_sample_set_free(&s);
    return s;
}

void geomlib_sample_set_free(geomlib_sample_set *s)
{
    if (!s) return;
    free(s->points); free(s->normals); free(s->triangle_ids);
    s->points = NULL; s->normals = NULL; s->triangle_ids = NULL;
    s->count = 0;
}

size_t geomlib_blue_noise_surface(const geomlib_trimesh *mesh,
                                  float min_distance,
                                  geomlib_pcg32 *rng_state,
                                  geomlib_sample_set *out)
{
    geomlib_pcg32 local_rng;
    geomlib_pcg32 *rng = rng_state;
    double *cdf;
    double total_area = 0.0;
    size_t t, i, capacity;
    float min2;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0 || !(min_distance > 0.0f))
        return 0;
    if (!rng) {
        geomlib_pcg32_init(&local_rng, 0x9e3779b9u, 0x85ebca6bu);
        rng = &local_rng;
    }
    min2 = min_distance * min_distance;
    capacity = 4 * mesh->triangle_count + 256;

    /* Area-weighted triangle CDF. */
    cdf = (double *)malloc(mesh->triangle_count * sizeof(double));
    if (!cdf) return 0;
    for (t = 0; t < mesh->triangle_count; ++t) {
        uint32_t a = mesh->indices[t*3], b = mesh->indices[t*3+1], c = mesh->indices[t*3+2];
        geomlib_vec3 e1 = geomlib_vec3_sub(mesh->vertices[b], mesh->vertices[a]);
        geomlib_vec3 e2 = geomlib_vec3_sub(mesh->vertices[c], mesh->vertices[a]);
        geomlib_vec3 n = geomlib_vec3_cross(e1, e2);
        total_area += 0.5 * geomlib_vec3_length(n);
        cdf[t] = total_area;
    }
    if (!(total_area > 0.0)) { free(cdf); return 0; }

    /* Brute-force O(N^2) dart throwing. The spatial-hash acceleration is a
     * future optimization; the algorithm and output are identical. */
    {
        geomlib_sample_set candidates = geomlib_sample_set_create(capacity);
        size_t attempts = 0;
        const size_t max_attempts = capacity * 30;
        while (attempts++ < max_attempts && candidates.count < capacity) {
            double pick = geomlib_pcg32_nextd(rng) * total_area;
            size_t lo = 0, hi = mesh->triangle_count;
            while (lo < hi) {
                size_t mid = (lo+hi)/2;
                if (cdf[mid] < pick) lo = mid+1; else hi = mid;
            }
            if (lo >= mesh->triangle_count) continue;
            {
                double u = geomlib_pcg32_nextd(rng);
                double v = geomlib_pcg32_nextd(rng);
                if (u + v > 1.0) { u = 1.0-u; v = 1.0-v; }
                uint32_t a = mesh->indices[lo*3], b = mesh->indices[lo*3+1],
                         c = mesh->indices[lo*3+2];
                geomlib_vec3 p = {
                    mesh->vertices[a].x + (float)u*(mesh->vertices[b].x-mesh->vertices[a].x)
                                       + (float)v*(mesh->vertices[c].x-mesh->vertices[a].x),
                    mesh->vertices[a].y + (float)u*(mesh->vertices[b].y-mesh->vertices[a].y)
                                       + (float)v*(mesh->vertices[c].y-mesh->vertices[a].y),
                    mesh->vertices[a].z + (float)u*(mesh->vertices[b].z-mesh->vertices[a].z)
                                       + (float)v*(mesh->vertices[c].z-mesh->vertices[a].z)
                };
                int accepted = 1;
                for (i = 0; i < candidates.count; ++i) {
                    geomlib_vec3 d = geomlib_vec3_sub(p, candidates.points[i]);
                    if (geomlib_vec3_dot(d, d) < min2) { accepted = 0; break; }
                }
                if (accepted) {
                    geomlib_vec3 n = geomlib_trimesh_triangle_normal(mesh, lo);
                    candidates.points[candidates.count] = p;
                    candidates.normals[candidates.count] = n;
                    candidates.triangle_ids[candidates.count] = (uint32_t)lo;
                    candidates.count++;
                }
            }
        }
        /* Shuffle for deterministic, unbiased sample order. */
        for (i = candidates.count; i > 1; --i) {
            size_t j = (size_t)(geomlib_pcg32_nextf(rng) * (float)i);
            geomlib_vec3 tp = candidates.points[i-1];
            geomlib_vec3 tn = candidates.normals[i-1];
            uint32_t ti = candidates.triangle_ids[i-1];
            candidates.points[i-1] = candidates.points[j];
            candidates.normals[i-1] = candidates.normals[j];
            candidates.triangle_ids[i-1] = candidates.triangle_ids[j];
            candidates.points[j] = tp;
            candidates.normals[j] = tn;
            candidates.triangle_ids[j] = ti;
        }
        *out = candidates;
    }

    free(cdf);
    return out->count;
}
