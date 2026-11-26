/*
 * geomlib - Ambient occlusion (see geomlib_ambient_occlusion.h)
 *
 * Clean-room Fibonacci hemisphere sampling with BVH occlusion rays.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_ambient_occlusion.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int geomlib_mesh_ambient_occlusion(geomlib_trimesh *mesh, int samples,
                                   float max_distance, float *out)
{
    geomlib_bvh bvh;
    size_t i;
    int k;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->triangle_count || !mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;
    if (samples<1) samples=16;
    if (samples>2048) samples=2048;
    if (!(max_distance>0.0f)) max_distance=1e30f;

    if (geomlib_trimesh_compute_normals(mesh)!=GEOMLIB_MESH_SUCCESS ||
        !mesh->normals) return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_bvh_build(&bvh,mesh,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) return GEOMLIB_MESH_ERROR_INVALID;

    for (i=0;i<mesh->vertex_count;++i) {
        geomlib_vec3 n=geomlib_vec3_normalize(mesh->normals[i]);
        geomlib_vec3 h,t1,t2;
        int blocked=0;
        if (geomlib_vec3_length(n)<0.5f) n=geomlib_vec3_create(0,1,0);
        h=fabsf(n.y)<0.9f ? geomlib_vec3_create(0,1,0)
                          : geomlib_vec3_create(1,0,0);
        t1=geomlib_vec3_normalize(geomlib_vec3_cross(h,n));
        t2=geomlib_vec3_cross(n,t1);
        for (k=0;k<samples;++k) {
            double z=(k+0.5)/(double)samples;
            double r=sqrt(fmax(0.0,1.0-z*z));
            double phi=k*2.399963229728653;
            geomlib_vec3 dir=geomlib_vec3_add(
                geomlib_vec3_mul(t1,(float)(r*cos(phi))),
                geomlib_vec3_add(geomlib_vec3_mul(t2,(float)(r*sin(phi))),
                                 geomlib_vec3_mul(n,(float)z)));
            geomlib_ray ray=geomlib_ray_create(
                geomlib_vec3_add(mesh->vertices[i],geomlib_vec3_mul(n,1e-4f)),
                dir);
            ray.t_min=1e-4f;
            ray.t_max=max_distance;
            if (geomlib_bvh_test_occlusion(&bvh,&ray)) blocked++;
        }
        out[i]=1.0f-(float)blocked/samples;
    }
    geomlib_bvh_free(&bvh);
    return GEOMLIB_MESH_SUCCESS;
}
