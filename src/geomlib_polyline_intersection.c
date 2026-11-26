/*
 * geomlib - Polyline intersection (see geomlib_polyline_intersection.h)
 *
 * Clean-room BVH segment tracing. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polyline_intersection.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_polyline_hits geomlib_polyline_hits_create(void)
{
    geomlib_polyline_hits h;
    memset(&h, 0, sizeof(h));
    return h;
}

void geomlib_polyline_hits_free(geomlib_polyline_hits *h)
{
    if (!h) return;
    free(h->points);
    memset(h, 0, sizeof(*h));
}

static int add_point(geomlib_polyline_hits *h, geomlib_vec3 p, size_t *cap)
{
    size_t k;
    for (k=0;k<h->count;++k)
        if (geomlib_vec3_length(geomlib_vec3_sub(h->points[k],p))<1e-5f) return 0;
    if (h->count==*cap) {
        size_t nc=*cap?(*cap)*2:32;
        geomlib_vec3 *np=(geomlib_vec3*)realloc(h->points,nc*sizeof(geomlib_vec3));
        if (!np) return -1;
        h->points=np; *cap=nc;
    }
    h->points[h->count++]=p;
    return 0;
}

int geomlib_polyline_mesh_intersection(const geomlib_vec3 *polyline, size_t n,
                                       const geomlib_trimesh *mesh,
                                       geomlib_polyline_hits *out)
{
    geomlib_polyline_hits r=geomlib_polyline_hits_create();
    geomlib_bvh bvh;
    size_t cap=0, i;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !polyline || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->triangle_count || n<2) return rc;
    if (geomlib_bvh_build(&bvh,mesh,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) return rc;

    for (i=0;i+1<n;++i) {
        geomlib_vec3 p=polyline[i],q=polyline[i+1];
        geomlib_vec3 delta=geomlib_vec3_sub(q,p);
        float len=geomlib_vec3_length(delta);
        geomlib_vec3 dir,o=p;
        int iter;
        if (!(len>1e-8f)) continue;
        dir=geomlib_vec3_mul(delta,1.0f/len);
        for (iter=0;iter<64;++iter) {
            geomlib_ray ray=geomlib_ray_create(o,dir);
            geomlib_ray_hit hit;
            ray.t_min=1e-5f;
            ray.t_max=len-1e-6f;
            if (ray.t_max<=ray.t_min ||
                !geomlib_bvh_trace_ray(&bvh,&ray,&hit)) break;
            if (add_point(&r,hit.position,&cap)<0) goto fail;
            o=geomlib_vec3_add(hit.position,geomlib_vec3_mul(dir,1e-4f));
        }
    }
    geomlib_bvh_free(&bvh);
    *out=r;
    return GEOMLIB_MESH_SUCCESS;

fail:
    geomlib_bvh_free(&bvh);
    geomlib_polyline_hits_free(&r);
    return rc;
}
