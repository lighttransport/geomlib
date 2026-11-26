/*
 * geomlib - Geodesic diameter (see geomlib_geodesic_diameter.h)
 *
 * Clean-room exact and double-sweep variants using the public Dijkstra API.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_diameter.h"
#include "geomlib/geomlib_geodesic_dijkstra.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int farthest_from(const geomlib_trimesh *m, uint32_t src,
                         uint32_t *out_v, float *out_d)
{
    geomlib_geodesic_distances d=geomlib_geodesic_create(m->vertex_count);
    size_t i;
    uint32_t best=src;
    float bd=-1.0f;
    int rc=geomlib_geodesic_dijkstra_from_vertex(m,src,&d);
    if (rc!=GEOMLIB_MESH_SUCCESS) { geomlib_geodesic_free(&d); return rc; }
    for (i=0;i<m->vertex_count;++i)
        if (d.distances[i]>bd) { bd=d.distances[i]; best=(uint32_t)i; }
    geomlib_geodesic_free(&d);
    *out_v=best; *out_d=bd;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_geodesic_diameter_exact(const geomlib_trimesh *mesh,
                                    geomlib_geodesic_diameter *out)
{
    geomlib_geodesic_distances d;
    size_t s,i;
    int rc=GEOMLIB_MESH_ERROR_INVALID;
    if (!out||!mesh||!mesh->vertices||!mesh->triangle_count) return rc;
    d=geomlib_geodesic_create(mesh->vertex_count);
    out->source=0; out->target=0; out->distance=0.0f;
    for (s=0;s<mesh->vertex_count;++s) {
        rc=geomlib_geodesic_dijkstra_from_vertex(mesh,(uint32_t)s,&d);
        if (rc!=GEOMLIB_MESH_SUCCESS) { geomlib_geodesic_free(&d); return rc; }
        for (i=0;i<mesh->vertex_count;++i)
            if (d.distances[i]>out->distance) {
                out->distance=d.distances[i];
                out->source=(uint32_t)s; out->target=(uint32_t)i;
            }
    }
    geomlib_geodesic_free(&d);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_geodesic_diameter_sweep(const geomlib_trimesh *mesh,
                                    geomlib_geodesic_diameter *out)
{
    uint32_t a=0,b;
    float d1,d2;
    int rc;
    if (!out||!mesh||!mesh->vertices||!mesh->triangle_count)
        return GEOMLIB_MESH_ERROR_INVALID;
    rc=farthest_from(mesh,0,&a,&d1);
    if (rc!=GEOMLIB_MESH_SUCCESS) return rc;
    rc=farthest_from(mesh,a,&b,&d2);
    if (rc!=GEOMLIB_MESH_SUCCESS) return rc;
    out->source=a; out->target=b; out->distance=d2;
    return GEOMLIB_MESH_SUCCESS;
}
