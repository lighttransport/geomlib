/*
 * geomlib - Geodesic ball (see geomlib_geodesic_ball.h)
 *
 * Clean-room Dijkstra radius query on the triangle-edge graph.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_geodesic_ball.h"
#include "geomlib/geomlib_geodesic_dijkstra.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_geodesic_ball geomlib_geodesic_ball_create(void)
{
    geomlib_geodesic_ball b;
    memset(&b,0,sizeof(b));
    return b;
}

void geomlib_geodesic_ball_free(geomlib_geodesic_ball *b)
{
    if (!b) return;
    free(b->mask);
    free(b->vertices);
    memset(b,0,sizeof(*b));
}

int geomlib_extract_geodesic_ball(const geomlib_trimesh *mesh,
                                  uint32_t source, float radius,
                                  geomlib_geodesic_ball *out)
{
    geomlib_geodesic_ball r=geomlib_geodesic_ball_create();
    geomlib_geodesic_distances d;
    size_t i;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->vertex_count || source>=mesh->vertex_count ||
        !(radius>=0.0f)) return rc;
    d=geomlib_geodesic_create(mesh->vertex_count);
    rc=geomlib_geodesic_dijkstra_from_vertex(mesh,source,&d);
    if (rc!=GEOMLIB_MESH_SUCCESS) { geomlib_geodesic_free(&d); return rc; }

    r.mask=(uint8_t*)calloc(mesh->vertex_count,1);
    r.vertices=(uint32_t*)malloc(mesh->vertex_count*sizeof(uint32_t));
    if (!r.mask||!r.vertices) {
        geomlib_geodesic_free(&d);
        geomlib_geodesic_ball_free(&r);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i=0;i<mesh->vertex_count;++i)
        if (d.distances[i]<=radius+1e-6f) {
            r.mask[i]=1;
            r.vertices[r.count++]=(uint32_t)i;
        }
    geomlib_geodesic_free(&d);
    *out=r;
    return GEOMLIB_MESH_SUCCESS;
}
