/*
 * geomlib - Surface Voronoi segmentation (see geomlib_surface_voronoi.h)
 *
 * Clean-room implementation: one Dijkstra per site, assign each vertex to the
 * strictly nearest site (ties keep the lower site index).
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_surface_voronoi.h"
#include "geomlib/geomlib_geodesic_dijkstra.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_surface_voronoi geomlib_surface_voronoi_create(void)
{
    geomlib_surface_voronoi v;
    memset(&v, 0, sizeof(v));
    return v;
}

void geomlib_surface_voronoi_free(geomlib_surface_voronoi *v)
{
    if (!v) return;
    free(v->labels);
    free(v->sizes);
    memset(v, 0, sizeof(*v));
}

int geomlib_surface_voronoi_compute(const geomlib_trimesh *mesh,
                                    const uint32_t *sites, size_t nsites,
                                    geomlib_surface_voronoi *out)
{
    geomlib_surface_voronoi result = geomlib_surface_voronoi_create();
    geomlib_geodesic_distances d;
    float *min_d=NULL;
    size_t i, k;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->vertex_count || !sites || !nsites) return rc;
    for (i=0;i<nsites;++i)
        if (sites[i]>=mesh->vertex_count) return rc;

    result.labels=(uint32_t*)malloc(mesh->vertex_count*sizeof(uint32_t));
    result.sizes=(size_t*)calloc(nsites,sizeof(size_t));
    min_d=(float*)malloc(mesh->vertex_count*sizeof(float));
    d=geomlib_geodesic_create(mesh->vertex_count);
    if (!result.labels || !result.sizes || !min_d || !d.distances) {
        geomlib_geodesic_free(&d);
        free(min_d);
        geomlib_surface_voronoi_free(&result);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i=0;i<mesh->vertex_count;++i) min_d[i]=FLT_MAX;

    for (k=0;k<nsites;++k) {
        rc=geomlib_geodesic_dijkstra(mesh,&sites[k],1,&d);
        if (rc!=GEOMLIB_MESH_SUCCESS) {
            geomlib_geodesic_free(&d);
            free(min_d);
            geomlib_surface_voronoi_free(&result);
            return rc;
        }
        for (i=0;i<mesh->vertex_count;++i)
            if (d.distances[i]<min_d[i]) {
                min_d[i]=d.distances[i];
                result.labels[i]=(uint32_t)k;
            }
    }
    for (i=0;i<mesh->vertex_count;++i) result.sizes[result.labels[i]]++;

    geomlib_geodesic_free(&d);
    free(min_d);
    result.vertex_count=mesh->vertex_count;
    result.site_count=nsites;
    *out=result;
    return GEOMLIB_MESH_SUCCESS;
}
