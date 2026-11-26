/*
 * geomlib - Approximate convex decomposition (see header)
 *
 * Clean-room combination of VSA region clustering and convex hull extraction.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_convex_decomposition.h"
#include "geomlib/geomlib_vsa.h"
#include "geomlib/geomlib_collision.h"

#include <stdlib.h>
#include <string.h>

geomlib_convex_decomposition geomlib_convex_decomposition_create(void)
{
    geomlib_convex_decomposition d;
    memset(&d, 0, sizeof(d));
    return d;
}

void geomlib_convex_decomposition_free(geomlib_convex_decomposition *d)
{
    size_t i;
    if (!d) return;
    for (i = 0; i < d->piece_count; ++i) geomlib_trimesh_free(&d->pieces[i]);
    free(d->pieces);
    memset(d, 0, sizeof(*d));
}

int geomlib_approx_convex_decomposition(const geomlib_trimesh *mesh,
                                        size_t region_count,
                                        geomlib_convex_decomposition *out)
{
    geomlib_vsa_result clustering = geomlib_vsa_result_create();
    geomlib_trimesh ignored = {0};
    geomlib_trimesh *pieces;
    size_t r, t, k;
    int rc;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (region_count < 1) region_count = 1;
    if (region_count > mesh->triangle_count) region_count = mesh->triangle_count;

    rc = geomlib_vsa_simplify(mesh, region_count, 5, &clustering, &ignored);
    if (rc != GEOMLIB_MESH_SUCCESS) {
        geomlib_vsa_result_free(&clustering);
        return rc;
    }
    geomlib_trimesh_free(&ignored);

    pieces = (geomlib_trimesh *)calloc(region_count, sizeof(geomlib_trimesh));
    if (!pieces) {
        geomlib_vsa_result_free(&clustering);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (r = 0; r < region_count; ++r) {
        geomlib_vec3 *points = (geomlib_vec3 *)malloc(mesh->triangle_count*3*sizeof(geomlib_vec3));
        size_t n = 0;
        if (!points) {
            for (k = 0; k < r; ++k) geomlib_trimesh_free(&pieces[k]);
            free(pieces);
            geomlib_vsa_result_free(&clustering);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (t = 0; t < mesh->triangle_count; ++t) {
            if (clustering.region_of_triangle[t] != r) continue;
            for (k = 0; k < 3; ++k) points[n++] = mesh->vertices[mesh->indices[t*3+k]];
        }
        if (n >= 4) {
            geomlib_trimesh temp = {0};
            size_t idx;
            if (geomlib_trimesh_init(&temp, n, 0) == GEOMLIB_MESH_SUCCESS) {
                for (idx = 0; idx < n; ++idx)
                    geomlib_trimesh_add_vertex(&temp, points[idx], NULL, NULL);
                geomlib_compute_convex_hull(&temp, &pieces[r]);
            }
            geomlib_trimesh_free(&temp);
        }
        free(points);
    }
    geomlib_vsa_result_free(&clustering);
    out->pieces = pieces;
    out->piece_count = region_count;
    return GEOMLIB_MESH_SUCCESS;
}
