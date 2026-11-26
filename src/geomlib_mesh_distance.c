/*
 * geomlib - Mesh-mesh minimum distance (see geomlib_mesh_distance.h)
 *
 * Clean-room dual BVH vertex-to-surface scan. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_mesh_distance.h"
#include "geomlib/geomlib_bvh.h"
#include "geomlib/geomlib_intersection.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_mesh_mesh_distance(const geomlib_trimesh *a,
                               const geomlib_trimesh *b,
                               geomlib_mesh_distance *out)
{
    geomlib_bvh bvha, bvhb;
    size_t i;
    double best=DBL_MAX;
    geomlib_vec3 best_a={0,0,0}, best_b={0,0,0};
    int have_a=0,have_b=0;

    if (!out || !a || !b || !a->vertices || !b->vertices ||
        !a->triangle_count || !b->triangle_count) return GEOMLIB_MESH_ERROR_INVALID;
    /* Intersecting surfaces have zero separation. */
    {
        geomlib_intersection inter=geomlib_mesh_mesh_intersection(a,b);
        if (inter.segment_count>0) {
            geomlib_intersection_free(&inter);
            out->point_a=geomlib_vec3_create(0,0,0);
            out->point_b=geomlib_vec3_create(0,0,0);
            out->distance=0.0f;
            return GEOMLIB_MESH_SUCCESS;
        }
        geomlib_intersection_free(&inter);
    }
    if (geomlib_bvh_build(&bvha,a,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_bvh_build(&bvhb,b,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) { geomlib_bvh_free(&bvha);
                                return GEOMLIB_MESH_ERROR_INVALID; }

    for (i=0;i<a->vertex_count;++i) {
        geomlib_vec3 cp; float d;
        if (!geomlib_bvh_closest_point(&bvhb,a->vertices[i],&cp,&d,NULL)) continue;
        if (d<best) { best=d; best_a=a->vertices[i]; best_b=cp; have_a=1; }
    }
    for (i=0;i<b->vertex_count;++i) {
        geomlib_vec3 cp; float d;
        if (!geomlib_bvh_closest_point(&bvha,b->vertices[i],&cp,&d,NULL)) continue;
        if (d<best) { best=d; best_a=cp; best_b=b->vertices[i]; have_b=1; }
    }
    (void)have_a; (void)have_b;
    geomlib_bvh_free(&bvha);
    geomlib_bvh_free(&bvhb);
    if (best==DBL_MAX) return GEOMLIB_MESH_ERROR_INVALID;
    out->point_a=best_a;
    out->point_b=best_b;
    out->distance=(float)best;
    return GEOMLIB_MESH_SUCCESS;
}
