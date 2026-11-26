/*
 * geomlib - 2D straight skeleton (see geomlib_straight_skeleton.h)
 *
 * Clean-room implementation of the wavefront/offset propagation method:
 * each polygon vertex moves inward along its angle bisector at the speed
 * required for unit-speed edge offsets. When an edge collapses, its vertices
 * merge. The traced vertex paths form the straight skeleton.
 *
 * Reference (paper only): Aichholzer et al., "A Straight Skeleton of a
 * Simple Polygon", 1995; Felkel & Obdrzalek, "Straight Skeleton
 * Implementation", 1998.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_straight_skeleton.h"

#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_straight_skeleton geomlib_straight_skeleton_create(void)
{
    geomlib_straight_skeleton s;
    memset(&s, 0, sizeof(s));
    return s;
}

void geomlib_straight_skeleton_free(geomlib_straight_skeleton *s)
{
    if (!s) return;
    free(s->points);
    free(s->edges);
    free(s->times);
    memset(s, 0, sizeof(*s));
}

typedef struct { geomlib_vec2 p; float t; } sample;

static float vec2_len(geomlib_vec2 v)
{
    return sqrtf(v.x*v.x + v.y*v.y);
}

int geomlib_compute_straight_skeleton(const geomlib_polygon2d *polygon,
                                      geomlib_straight_skeleton *out)
{
    size_t n, i, j, nsamples = 0, step;
    geomlib_vec2 *front, *next_front;
    sample **paths;
    size_t *path_lengths, *path_caps;
    bool *alive;
    float time = 0.0f;
    const float dt = 0.01f;
    const size_t max_steps = 600;
    float area0;

    if (!polygon || !out || !polygon->vertices || polygon->vertex_count < 3)
        return GEOMLIB_MESH_ERROR_INVALID;
    n = polygon->vertex_count;
    area0 = geomlib_polygon2d_area(polygon);
    if (!(fabsf(area0) > 1e-9f)) return GEOMLIB_MESH_ERROR_INVALID;

    front = (geomlib_vec2 *)malloc(n*sizeof(geomlib_vec2));
    next_front = (geomlib_vec2 *)malloc(n*sizeof(geomlib_vec2));
    paths = (sample **)calloc(n, sizeof(sample*));
    path_lengths = (size_t *)calloc(n, sizeof(size_t));
    path_caps = (size_t *)malloc(n*sizeof(size_t));
    alive = (bool *)calloc(n, sizeof(bool));
    if (!front || !next_front || !paths || !path_lengths || !path_caps || !alive) {
        free(front); free(next_front); free(paths);
        free(path_lengths); free(path_caps); free(alive);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < n; ++i) {
        front[i] = polygon->vertices[i];
        alive[i] = true;
        path_caps[i] = 32;
        paths[i] = (sample *)malloc(path_caps[i]*sizeof(sample));
        if (!paths[i]) {
            for (j = 0; j < i; ++j) free(paths[j]);
            free(front); free(next_front); free(paths);
            free(path_lengths); free(path_caps); free(alive);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        paths[i][path_lengths[i]++] = (sample){front[i], 0.0f};
    }

    for (step = 0; step < max_steps; ++step) {
        int any_alive = 0;
        float area = 0.0f;
        for (i = 0; i < n; ++i) {
            if (!alive[i]) { next_front[i] = front[i]; continue; }
            any_alive = 1;
            geomlib_vec2 prev = front[(i+n-1)%n];
            geomlib_vec2 cur = front[i];
            geomlib_vec2 next = front[(i+1)%n];
            geomlib_vec2 e1 = geomlib_vec2_sub(cur, prev);
            geomlib_vec2 e2 = geomlib_vec2_sub(next, cur);
            float l1 = vec2_len(e1), l2 = vec2_len(e2);
            /* Inward normal for CCW polygon: (-dy, dx). */
            geomlib_vec2 n1 = {-e1.y/(l1>1e-9f?l1:1.0f), e1.x/(l1>1e-9f?l1:1.0f)};
            geomlib_vec2 n2 = {-e2.y/(l2>1e-9f?l2:1.0f), e2.x/(l2>1e-9f?l2:1.0f)};
            geomlib_vec2 bis = geomlib_vec2_add(n1, n2);
            float bl = vec2_len(bis);
            if (bl < 1e-6f) { next_front[i] = cur; continue; }
            bis.x /= bl; bis.y /= bl;
            /* Offset both incident lines by dt; solve d·v = dt. */
            {
                float denom = bis.x*n1.x + bis.y*n1.y;
                if (fabsf(denom) < 1e-4f) { next_front[i] = cur; continue; }
                float s = dt / denom;
                next_front[i] = geomlib_vec2_add(cur, geomlib_vec2_mul(bis, s));
            }
        }
        if (!any_alive) break;
        memcpy(front, next_front, n*sizeof(geomlib_vec2));

        /* Collapse edges that have become very short or flipped. */
        for (i = 0; i < n; ++i) {
            if (!alive[i] || !alive[(i+1)%n]) continue;
            geomlib_vec2 d = geomlib_vec2_sub(front[(i+1)%n], front[i]);
            if (vec2_len(d) > 2.0f*dt) continue;
            /* Merge the two vertices: the survivor keeps its path. */
            geomlib_vec2 merged = geomlib_vec2_mul(
                geomlib_vec2_add(front[i], front[(i+1)%n]), 0.5f);
            size_t survivor = i;
            size_t victim = (i+1)%n;
            front[survivor] = merged;
            alive[victim] = false;
            if (path_lengths[survivor] == path_caps[survivor]) {
                size_t nc = path_caps[survivor]*2;
                sample *np = (sample *)realloc(paths[survivor], nc*sizeof(sample));
                if (np) { paths[survivor] = np; path_caps[survivor] = nc; }
            }
            if (path_lengths[survivor] < path_caps[survivor])
                paths[survivor][path_lengths[survivor]++] = (sample){merged, time + dt};
        }

        time += dt;
        /* Area of the (possibly partially collapsed) front. */
        {
            geomlib_polygon2d poly = {front, n, n};
            area = fabsf(geomlib_polygon2d_area(&poly));
        }
        if (area < 1e-4f) break;

        /* Record traced paths for surviving vertices. */
        for (i = 0; i < n; ++i) {
            if (!alive[i]) continue;
            if (path_lengths[i] == path_caps[i]) {
                size_t nc = path_caps[i]*2;
                sample *np = (sample *)realloc(paths[i], nc*sizeof(sample));
                if (np) { paths[i] = np; path_caps[i] = nc; }
            }
            if (path_lengths[i] < path_caps[i])
                paths[i][path_lengths[i]++] = (sample){front[i], time};
        }
    }

    /* Flatten the traced paths into the skeleton graph. */
    for (i = 0; i < n; ++i) nsamples += path_lengths[i];
    if (nsamples == 0) {
        for (i = 0; i < n; ++i) free(paths[i]);
        free(front); free(next_front); free(paths);
        free(path_lengths); free(path_caps); free(alive);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    {
        geomlib_straight_skeleton result = geomlib_straight_skeleton_create();
        size_t idx = 0, ecount = 0;
        result.points = (geomlib_vec2 *)malloc(nsamples*sizeof(geomlib_vec2));
        result.times = (float *)malloc(nsamples*sizeof(float));
        result.edges = (uint32_t *)malloc(nsamples*2*sizeof(uint32_t));
        if (!result.points || !result.times || !result.edges) {
            geomlib_straight_skeleton_free(&result);
        } else {
            for (i = 0; i < n; ++i) {
                size_t base = idx;
                for (j = 0; j < path_lengths[i]; ++j) {
                    result.points[idx] = paths[i][j].p;
                    result.times[idx] = paths[i][j].t;
                    if (j > 0) {
                        result.edges[ecount*2] = (uint32_t)(base+j-1);
                        result.edges[ecount*2+1] = (uint32_t)(base+j);
                        ecount++;
                    }
                    idx++;
                }
            }
            result.point_count = idx;
            result.edge_count = ecount;
            *out = result;
        }
    }

    for (i = 0; i < n; ++i) free(paths[i]);
    free(front); free(next_front); free(paths);
    free(path_lengths); free(path_caps); free(alive);
    return GEOMLIB_MESH_SUCCESS;
}
