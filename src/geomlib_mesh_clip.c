/*
 * geomlib - Plane clipping (see geomlib_mesh_clip.h)
 *
 * Clean-room per-triangle half-space clipping with cap reconstruction.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_mesh_clip.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int geomlib_clip_mesh_plane(const geomlib_trimesh *mesh,
                            geomlib_vec3 plane_point, geomlib_vec3 plane_normal,
                            int keep_positive, geomlib_trimesh *out)
{
    geomlib_vec3 n,*verts=NULL,*cut=NULL;
    double *d=NULL;
    size_t i,k,ncut=0,ncap=0;
    int sgn;

    if (!out || !mesh || !mesh->vertices || !mesh->indices) return GEOMLIB_MESH_ERROR_INVALID;
    memset(out,0,sizeof(*out));
    n=plane_normal;
    {
        float l=geomlib_vec3_length(n);
        if (!(l>1e-12f)) return GEOMLIB_MESH_ERROR_INVALID;
        n=geomlib_vec3_mul(n,1.0f/l);
    }
    sgn=keep_positive ? 1 : -1;
    d=(double*)malloc(mesh->vertex_count*sizeof(double));
    verts=(geomlib_vec3*)malloc(mesh->vertex_count*sizeof(geomlib_vec3));
    cut=(geomlib_vec3*)malloc(6*mesh->triangle_count*sizeof(geomlib_vec3));
    if (!d||!verts||!cut) { free(d);free(verts);free(cut); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i=0;i<mesh->vertex_count;++i) {
        verts[i]=mesh->vertices[i];
        d[i]=sgn*geomlib_vec3_dot(n,geomlib_vec3_sub(mesh->vertices[i],plane_point));
    }
    if (geomlib_trimesh_init(out,mesh->vertex_count+8,mesh->triangle_count+32)!=
        GEOMLIB_MESH_SUCCESS) { free(d);free(verts);free(cut); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i=0;i<mesh->vertex_count;++i)
        geomlib_trimesh_add_vertex(out,verts[i],NULL,NULL);

    for (i=0;i<mesh->triangle_count;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        geomlib_vec3 poly[4];
        size_t np=0;
        int inside[3];
        for (k=0;k<3;++k) inside[k]=d[v[k]]>=-1e-10;
        for (k=0;k<3;++k) {
            size_t nx=(k+1)%3;
            if (inside[k]) poly[np++]=verts[v[k]];
            if (inside[k]!=inside[nx]) {
                double t=d[v[k]]/(d[v[k]]-d[v[nx]]);
                geomlib_vec3 p=geomlib_vec3_add(verts[v[k]],
                    geomlib_vec3_mul(geomlib_vec3_sub(verts[v[nx]],verts[v[k]]),(float)t));
                poly[np++]=p;
                if (ncut<6*mesh->triangle_count) cut[ncut++]=p;
            }
        }
        if (np<3) continue;
        {
            size_t base=out->vertex_count;
            for (k=0;k<np;++k) geomlib_trimesh_add_vertex(out,poly[k],NULL,NULL);
            for (k=1;k+1<np;++k)
                geomlib_trimesh_add_triangle(out,
                    (uint32_t)base,(uint32_t)(base+k),(uint32_t)(base+k+1));
        }
    }

    /* Cap: deduplicate cut points, sort by angle in the plane, fan triangulate. */
    if (ncut>=3) {
        geomlib_vec3 center={0,0,0};
        geomlib_vec3 u,w;
        size_t m;
        for (i=0;i<ncut;++i) center=geomlib_vec3_add(center,cut[i]);
        center=geomlib_vec3_mul(center,1.0f/(float)ncut);
        {
            geomlib_vec3 h=fabsf(n.x)<0.9f?geomlib_vec3_create(1,0,0)
                                          :geomlib_vec3_create(0,1,0);
            u=geomlib_vec3_normalize(geomlib_vec3_cross(n,h));
            w=geomlib_vec3_cross(n,u);
        }
        {
            /* Simple insertion sort by angle. */
            for (i=1;i<ncut;++i) {
                geomlib_vec3 key=cut[i];
                double ak=atan2(geomlib_vec3_dot(geomlib_vec3_sub(key,center),w),
                                geomlib_vec3_dot(geomlib_vec3_sub(key,center),u));
                size_t j=i;
                while (j>0) {
                    geomlib_vec3 prev=cut[j-1];
                    double ap=atan2(geomlib_vec3_dot(geomlib_vec3_sub(prev,center),w),
                                    geomlib_vec3_dot(geomlib_vec3_sub(prev,center),u));
                    if (ap<=ak) break;
                    cut[j]=prev; --j;
                }
                cut[j]=key;
            }
        }
        /* Deduplicate consecutive coincident points. */
        m=0;
        for (i=0;i<ncut;++i) {
            if (m && geomlib_vec3_length(geomlib_vec3_sub(cut[m-1],cut[i]))<1e-5f) continue;
            cut[m++]=cut[i];
        }
        ncap=m;
        if (ncap>=3) {
            size_t base=out->vertex_count;
            for (i=0;i<ncap;++i) geomlib_trimesh_add_vertex(out,cut[i],NULL,NULL);
            for (i=1;i+1<ncap;++i) {
                size_t a=base,b=base+i,c=base+i+1;
                if (!keep_positive) { size_t t=b; b=c; c=t; }
                geomlib_trimesh_add_triangle(out,(uint32_t)a,(uint32_t)b,(uint32_t)c);
            }
        }
    }
    free(d); free(verts); free(cut);
    return GEOMLIB_MESH_SUCCESS;
}
