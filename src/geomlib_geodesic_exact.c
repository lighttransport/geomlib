/*
 * geomlib - Continuous geodesic path length (see geomlib_geodesic_exact.h)
 *
 * Clean-room implementation of face-strip unfolding: find a shortest edge
 * path, unfold each traversed triangle across the shared edge by solving a
 * two-circle intersection, and measure the straight line between endpoints
 * in the final unfolded plane.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_exact.h"
#include "geomlib/geomlib_topology.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

typedef struct { double x, y; } p2;

static double dist3(const geomlib_trimesh *m, uint32_t a, uint32_t b)
{
    geomlib_vec3 d = geomlib_vec3_sub(m->vertices[a], m->vertices[b]);
    return sqrt(geomlib_vec3_dot(d, d));
}

/* Find a triangle containing both vertices, preferring one also containing
 * a third hint vertex (so consecutive faces stay adjacent). */
static size_t find_face(const geomlib_trimesh *m, uint32_t a, uint32_t b,
                        uint32_t hint)
{
    size_t t;
    for (t = 0; t < m->triangle_count; ++t) {
        uint32_t v0=m->indices[t*3], v1=m->indices[t*3+1], v2=m->indices[t*3+2];
        int has_a = (v0==a||v1==a||v2==a);
        int has_b = (v0==b||v1==b||v2==b);
        int has_h = (v0==hint||v1==hint||v2==hint);
        if (has_a && has_b && has_h) return t;
    }
    for (t = 0; t < m->triangle_count; ++t) {
        uint32_t v0=m->indices[t*3], v1=m->indices[t*3+1], v2=m->indices[t*3+2];
        if ((v0==a||v1==a||v2==a) && (v0==b||v1==b||v2==b)) return t;
    }
    return (size_t)-1;
}

int geomlib_geodesic_length_unfolding(const geomlib_trimesh *mesh,
                                      uint32_t source,
                                      uint32_t target,
                                      double *length)
{
    geomlib_vertex_adjacency adj = {0};
    uint32_t *path;
    double *d;
    bool *done;
    size_t V, n, i, j;
    p2 *planes;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !length || !mesh->vertices || !mesh->indices ||
        source >= mesh->vertex_count || target >= mesh->vertex_count)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (source == target) { *length = 0.0; return GEOMLIB_MESH_SUCCESS; }

    d = (double *)malloc(V*sizeof(double));
    done = (bool *)calloc(V, sizeof(bool));
    path = (uint32_t *)malloc((V+1)*sizeof(uint32_t));
    planes = (p2 *)malloc((V+1)*sizeof(p2));
    if (!d || !done || !path || !planes) {
        free(d); free(done); free(path); free(planes);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < V; ++i) d[i] = INFINITY;
    d[source] = 0.0;
    adj = geomlib_build_vertex_adjacency(mesh);
    if (!adj.neighbor_counts) {
        free(d); free(done); free(path); free(planes);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Dijkstra edge path. */
    for (i = 0; i < V; ++i) {
        size_t best = V;
        double bestd = INFINITY;
        for (j = 0; j < V; ++j)
            if (!done[j] && d[j] < bestd) { bestd = d[j]; best = j; }
        if (best == V || !isfinite(bestd)) break;
        done[best] = true;
        for (j = 0; j < adj.neighbor_counts[best]; ++j) {
            uint32_t nb = adj.neighbors[best][j];
            double w = dist3(mesh, (uint32_t)best, nb);
            if (w < 1e-9) continue; /* duplicate vertices */
            if (d[best] + w < d[nb]) d[nb] = d[best] + w;
        }
    }
    if (!isfinite(d[target])) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }

    /* Reconstruct path target -> source. */
    n = 0;
    {
        uint32_t cur = target;
        path[n++] = cur;
        while (cur != source) {
            uint32_t prev = UINT32_MAX;
            for (j = 0; j < adj.neighbor_counts[cur]; ++j) {
                uint32_t nb = adj.neighbors[cur][j];
                double w = dist3(mesh, cur, nb);
                if (w < 1e-9) continue;
                if (fabs(d[nb] + w - d[cur]) < 1e-6) { prev = nb; break; }
            }
            if (prev == UINT32_MAX) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
            cur = prev;
            path[n++] = cur;
        }
    }
    /* Reverse to source -> target. */
    for (i = 0; i < n/2; ++i) {
        uint32_t t = path[i]; path[i] = path[n-1-i]; path[n-1-i] = t;
    }

    /* Unfold the strip. Place the first edge along the x axis. */
    planes[0].x = 0.0; planes[0].y = 0.0;
    {
        double len = dist3(mesh, path[0], path[1]);
        planes[1].x = len; planes[1].y = 0.0;
    }
    for (i = 0; i + 2 < n; ++i) {
        uint32_t a = path[i], b = path[i+1], c = path[i+2];
        size_t f = find_face(mesh, a, b, c);
        if (f == (size_t)-1) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
        (void)f;
        {
            p2 A = planes[i], B = planes[i+1];
            double ux = B.x-A.x, uy = B.y-A.y;
            double L = sqrt(ux*ux+uy*uy);
            if (L < 1e-9) { planes[i+2] = B; continue; }
            double ac = dist3(mesh, a, c), bc = dist3(mesh, b, c);
            double t = (ac*ac - bc*bc + L*L) / (2.0*L);
            double h2 = ac*ac - t*t;
            double h = sqrt(h2 > 0.0 ? h2 : 0.0);
            p2 C;
            C.x = A.x + t*ux/L - h*uy/L;
            C.y = A.y + t*uy/L + h*ux/L;
            planes[i+2] = C;
        }
    }
    {
        double dx = planes[n-1].x - planes[0].x;
        double dy = planes[n-1].y - planes[0].y;
        *length = sqrt(dx*dx + dy*dy);
    }
    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_vertex_adjacency_free(&adj);
    free(d); free(done); free(path); free(planes);
    return rc;
}
