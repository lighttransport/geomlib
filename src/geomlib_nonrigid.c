/*
 * geomlib - Non-rigid ICP (see geomlib_nonrigid.h)
 *
 * Clean-room implementation: alternate closest-point correspondence with a
 * Jacobi update of a per-vertex displacement field blended with mesh-neighbour
 * averages. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_nonrigid.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_nonrigid_result geomlib_nonrigid_result_create(void)
{
    geomlib_nonrigid_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_nonrigid_result_free(geomlib_nonrigid_result *r)
{
    if (!r) return;
    free(r->displacement);
    memset(r, 0, sizeof(*r));
}

typedef struct { size_t *offsets; uint32_t *neighbors; } adjacency;

static void adjacency_free(adjacency *a)
{
    free(a->offsets); free(a->neighbors);
    memset(a, 0, sizeof(*a));
}

static int build_adjacency(const geomlib_trimesh *m, adjacency *a)
{
    size_t n=m->vertex_count, i;
    size_t *count;
    memset(a, 0, sizeof(*a));
    if (!n || !m->indices) return GEOMLIB_MESH_ERROR_INVALID;
    count = (size_t *)calloc(n, sizeof(size_t));
    a->offsets = (size_t *)calloc(n+1, sizeof(size_t));
    if (!count || !a->offsets) { free(count); adjacency_free(a); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < m->triangle_count; ++i) {
        uint32_t v[3]={m->indices[i*3],m->indices[i*3+1],m->indices[i*3+2]};
        int k;
        if (v[0]>=n||v[1]>=n||v[2]>=n) continue;
        for (k=0;k<3;++k) {
            count[v[k]]++;
            count[v[(k+1)%3]]++;
        }
    }
    for (i = 0; i < n; ++i) a->offsets[i+1]=a->offsets[i]+count[i];
    a->neighbors = (uint32_t *)malloc((a->offsets[n]?a->offsets[n]:1)*sizeof(uint32_t));
    if (!a->neighbors) { free(count); adjacency_free(a); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < n; ++i) count[i]=0;
    for (i = 0; i < m->triangle_count; ++i) {
        uint32_t v[3]={m->indices[i*3],m->indices[i*3+1],m->indices[i*3+2]};
        int k;
        if (v[0]>=n||v[1]>=n||v[2]>=n) continue;
        for (k=0;k<3;++k) {
            uint32_t x=v[k], y=v[(k+1)%3];
            a->neighbors[a->offsets[x]+count[x]++]=y;
            a->neighbors[a->offsets[y]+count[y]++]=x;
        }
    }
    free(count);
    return GEOMLIB_MESH_SUCCESS;
}

static float correspondence_rms(const geomlib_trimesh *src,
                                const geomlib_bvh *bvh,
                                const geomlib_vec3 *d)
{
    double sum=0.0;
    size_t i;
    for (i = 0; i < src->vertex_count; ++i) {
        geomlib_vec3 q=geomlib_vec3_add(src->vertices[i],d[i]);
        geomlib_vec3 cp;
        float dist;
        if (!geomlib_bvh_closest_point(bvh,q,&cp,&dist,NULL)) dist=FLT_MAX;
        if (!isfinite(dist)) dist=FLT_MAX;
        if (dist < FLT_MAX) sum += (double)dist*dist;
    }
    return src->vertex_count ? (float)sqrt(sum/src->vertex_count) : 0.0f;
}

int geomlib_nonrigid_icp(const geomlib_trimesh *source,
                         const geomlib_trimesh *target,
                         uint32_t max_iterations, float smoothness,
                         geomlib_nonrigid_result *out)
{
    geomlib_nonrigid_result r=geomlib_nonrigid_result_create();
    geomlib_bvh bvh;
    adjacency adj;
    geomlib_vec3 *d=NULL, *dn=NULL;
    size_t n, i, iter;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !source || !target || !source->vertices || !source->indices ||
        !target->vertices || !target->indices || !source->vertex_count ||
        !target->triangle_count) return rc;
    if (!(smoothness>=0.0f && smoothness<=1.0f)) smoothness=0.85f;
    if (!max_iterations) max_iterations=1;
    n=source->vertex_count;

    d=(geomlib_vec3*)calloc(n,sizeof(geomlib_vec3));
    dn=(geomlib_vec3*)calloc(n,sizeof(geomlib_vec3));
    if (!d||!dn) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
    rc=build_adjacency(source,&adj);
    if (rc!=GEOMLIB_MESH_SUCCESS) goto cleanup;
    if (geomlib_bvh_build(&bvh,target,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) { rc=GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }

    r.initial_rms=correspondence_rms(source,&bvh,d);
    for (iter=0;iter<max_iterations;++iter) {
        for (i=0;i<n;++i) {
            geomlib_vec3 q=geomlib_vec3_add(source->vertices[i],d[i]);
            geomlib_vec3 cp, goal=d[i];
            float dist;
            if (geomlib_bvh_closest_point(&bvh,q,&cp,&dist,NULL) &&
                isfinite(dist)) goal=geomlib_vec3_add(d[i],geomlib_vec3_sub(cp,q));
            {
                geomlib_vec3 avg=geomlib_vec3_create(0,0,0);
                size_t k,a=adj.offsets[i],b=adj.offsets[i+1];
                for (k=a;k<b;++k) avg=geomlib_vec3_add(avg,d[adj.neighbors[k]]);
                if (b>a) avg=geomlib_vec3_mul(avg,1.0f/(float)(b-a));
                dn[i]=geomlib_vec3_add(
                    geomlib_vec3_mul(goal,smoothness),
                    geomlib_vec3_mul(avg,1.0f-smoothness));
            }
        }
        { geomlib_vec3 *tmp=d; d=dn; dn=tmp; }
    }
    r.final_rms=correspondence_rms(source,&bvh,d);
    r.displacement=d; d=NULL;
    r.vertex_count=n;
    r.iterations=max_iterations;
    *out=r;
    rc=GEOMLIB_MESH_SUCCESS;

    geomlib_bvh_free(&bvh);
cleanup:
    adjacency_free(&adj);
    free(d); free(dn);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_nonrigid_result_free(&r);
    return rc;
}

int geomlib_nonrigid_apply(const geomlib_trimesh *source,
                           const geomlib_nonrigid_result *fit,
                           geomlib_trimesh *out)
{
    size_t i;
    if (!out || !source || !fit || !source->vertices ||
        fit->vertex_count != source->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_trimesh_init(out, source->vertex_count, source->triangle_count)!=
        GEOMLIB_MESH_SUCCESS) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i=0;i<source->vertex_count;++i) {
        geomlib_vec3 p=geomlib_vec3_add(source->vertices[i],fit->displacement[i]);
        if (geomlib_trimesh_add_vertex(out,p,NULL,NULL)<0) {
            geomlib_trimesh_free(out);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    if (source->indices)
        for (i=0;i<source->triangle_count;++i)
            if (geomlib_trimesh_add_triangle(out,
                    source->indices[i*3], source->indices[i*3+1],
                    source->indices[i*3+2])<0) {
                geomlib_trimesh_free(out);
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
    return GEOMLIB_MESH_SUCCESS;
}
