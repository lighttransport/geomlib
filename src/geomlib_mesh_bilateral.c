/*
 * geomlib - Bilateral mesh smoothing (see geomlib_mesh_bilateral.h)
 *
 * Clean-room implementation: k-d-tree neighbours, Gaussian spatial and normal
 * weights, simultaneous (Jacobi) position updates. C11, wasm32 safe.
 */

#include "geomlib/geomlib_mesh_bilateral.h"
#include "geomlib/geomlib_kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_mesh_bilateral_params geomlib_mesh_bilateral_params_default(void)
{
    geomlib_mesh_bilateral_params p;
    p.sigma_spatial=0.0f;
    p.sigma_normal=0.35f;
    p.neighbors=12;
    p.iterations=2;
    return p;
}

int geomlib_mesh_bilateral_smooth(geomlib_trimesh *mesh,
                                  const geomlib_mesh_bilateral_params *params)
{
    geomlib_mesh_bilateral_params p = params ? *params
                                             : geomlib_mesh_bilateral_params_default();
    geomlib_kdtree *tree=NULL;
    uint32_t *idx=NULL;
    float *dst=NULL, *normalized=NULL;
    geomlib_vec3 *next=NULL;
    size_t n, i, k, iter;
    float sigma;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !mesh->vertices || !mesh->vertex_count ||
        !mesh->indices || !mesh->triangle_count) return rc;
    n=mesh->vertex_count;
    if (p.neighbors<3) p.neighbors=3;
    if (p.neighbors > 64) p.neighbors=64;
    if (p.iterations<1) p.iterations=1;
    if (p.iterations>32) p.iterations=32;
    if (!(p.sigma_normal>0.0f)) p.sigma_normal=0.35f;

    if (geomlib_trimesh_compute_normals(mesh) != GEOMLIB_MESH_SUCCESS)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!mesh->normals) return GEOMLIB_MESH_ERROR_INVALID;
    normalized=(float*)malloc(n*3*sizeof(float));
    if (!normalized) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i=0;i<n;++i) {
        geomlib_vec3 v=mesh->normals[i];
        float l=geomlib_vec3_length(v);
        if (l>1e-12f) v=geomlib_vec3_mul(v,1.0f/l);
        else v=geomlib_vec3_create(0,0,1);
        normalized[i*3]=v.x; normalized[i*3+1]=v.y; normalized[i*3+2]=v.z;
    }
    {
        float mean=0.0f;
        for (i=0;i<mesh->triangle_count;++i) {
            uint32_t a=mesh->indices[i*3],b=mesh->indices[i*3+1],c=mesh->indices[i*3+2];
            if (a>=n||b>=n||c>=n) continue;
            mean+=geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[a],mesh->vertices[b]));
            mean+=geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[b],mesh->vertices[c]));
        }
        sigma = p.sigma_spatial>0.0f ? p.sigma_spatial
              : (mesh->triangle_count ? mean/(2.0f*mesh->triangle_count) : 1.0f);
        if (!(sigma>1e-9f)) sigma=1e-3f;
    }

    k=(size_t)p.neighbors;
    idx=(uint32_t*)malloc(k*sizeof(uint32_t));
    dst=(float*)malloc(k*sizeof(float));
    next=(geomlib_vec3*)malloc(n*sizeof(geomlib_vec3));
    if (!idx||!dst||!next) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }

    for (iter=0;iter<(size_t)p.iterations;++iter) {
        tree=geomlib_kdtree_build(mesh->vertices,n);
        if (!tree) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
        for (i=0;i<n;++i) {
            size_t j, found;
            double wsum=0.0, x=0.0,y=0.0,z=0.0;
            found=geomlib_kdtree_knn(tree,mesh->vertices[i],k,idx,dst);
            if (!found) { next[i]=mesh->vertices[i]; continue; }
            for (j=0;j<found;++j) {
                uint32_t v=idx[j];
                double ds=(double)dst[j]/sigma;
                double nx=(double)normalized[i*3]*normalized[v*3]
                         + (double)normalized[i*3+1]*normalized[v*3+1]
                         + (double)normalized[i*3+2]*normalized[v*3+2];
                double dn=(1.0-nx)/(double)(p.sigma_normal*p.sigma_normal);
                double w=exp(-ds*ds-dn);
                if (!isfinite(w)) continue;
                x+=w*mesh->vertices[v].x;
                y+=w*mesh->vertices[v].y;
                z+=w*mesh->vertices[v].z;
                wsum+=w;
            }
            next[i] = wsum>1e-18 ? geomlib_vec3_create((float)(x/wsum),
                                                        (float)(y/wsum),
                                                        (float)(z/wsum))
                                 : mesh->vertices[i];
        }
        geomlib_kdtree_free(tree); tree=NULL;
        memcpy(mesh->vertices,next,n*sizeof(geomlib_vec3));
    }
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_kdtree_free(tree);
    free(idx); free(dst); free(next); free(normalized);
    return rc;
}
