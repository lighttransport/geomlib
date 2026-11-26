/*
 * geomlib - Point-to-plane ICP (see geomlib_icp_plane.h)
 *
 * Clean-room implementation of the Chen-Medioni linearized least-squares
 * rigid fit. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_icp_plane.h"
#include "geomlib/geomlib_bvh.h"
#include "geomlib/geomlib_solver.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_icp_plane_params geomlib_icp_plane_params_default(void)
{
    geomlib_icp_plane_params p;
    p.max_iterations=30;
    p.tolerance=1e-4f;
    return p;
}

static geomlib_vec3 transform_point(geomlib_mat4 t, geomlib_vec3 p)
{
    geomlib_vec4 v=geomlib_mat4_mul_vec4(t,
        geomlib_vec4_create(p.x,p.y,p.z,1.0f));
    return geomlib_vec3_create(v.x,v.y,v.z);
}

static geomlib_mat4 rodrigues(geomlib_vec3 w, geomlib_vec3 t)
{
    geomlib_mat4 m=geomlib_mat4_identity();
    float angle=geomlib_vec3_length(w);
    float c,s;
    geomlib_vec3 axis;
    if (angle < 1e-9f) {
        c=1.0f; s=0.0f; axis=geomlib_vec3_create(0,0,1);
    } else {
        c=cosf(angle); s=sinf(angle); axis=geomlib_vec3_mul(w,1.0f/angle);
    }
    {
        float x=axis.x,y=axis.y,z=axis.z;
        float ic=1.0f-c;
        /* Column-major: m[col*4+row]. */
        m.m[0]=c+x*x*ic;     m.m[4]=x*y*ic-z*s;   m.m[8] =x*z*ic+y*s;
        m.m[1]=y*x*ic+z*s;   m.m[5]=c+y*y*ic;     m.m[9] =y*z*ic-x*s;
        m.m[2]=z*x*ic-y*s;   m.m[6]=z*y*ic+x*s;   m.m[10]=c+z*z*ic;
    }
    m.m[12]=t.x; m.m[13]=t.y; m.m[14]=t.z;
    return m;
}

int geomlib_icp_point_to_plane(const geomlib_trimesh *source,
                               const geomlib_trimesh *target,
                               const geomlib_icp_plane_params *params,
                               geomlib_mat4 *transform, float *final_rms)
{
    geomlib_icp_plane_params p=params ? *params : geomlib_icp_plane_params_default();
    geomlib_bvh bvh;
    geomlib_mat4 T=geomlib_mat4_identity();
    double prev=DBL_MAX;
    int iter, rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!source || !target || !transform || !source->vertices ||
        !target->vertices || !target->indices || !source->vertex_count ||
        !target->triangle_count) return rc;
    if (p.max_iterations<1) p.max_iterations=1;
    if (!(p.tolerance>0.0f)) p.tolerance=1e-4f;
    if (geomlib_bvh_build(&bvh,target,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) return rc;

    for (iter=0;iter<p.max_iterations;++iter) {
        double ata[36]={0}, atb[6]={0};
        double sum=0.0;
        size_t i;
        int solve_rc;
        for (i=0;i<source->vertex_count;++i) {
            geomlib_vec3 q0=transform_point(T,source->vertices[i]);
            geomlib_vec3 cp;
            float dist;
            uint32_t fid=0;
            geomlib_vec3 n=geomlib_vec3_create(0,0,1);
            double a[6], b;
            int k,l;
            if (!geomlib_bvh_closest_point(&bvh,q0,&cp,&dist,&fid)) continue;
            if (fid<target->triangle_count) {
                uint32_t i0=target->indices[fid*3],i1=target->indices[fid*3+1],
                         i2=target->indices[fid*3+2];
                if (i0<target->vertex_count&&i1<target->vertex_count&&i2<target->vertex_count)
                    n=geomlib_vec3_normalize(geomlib_vec3_cross(
                        geomlib_vec3_sub(target->vertices[i1],target->vertices[i0]),
                        geomlib_vec3_sub(target->vertices[i2],target->vertices[i0])));
            }
            {
                geomlib_vec3 r=geomlib_vec3_cross(q0,n);
                a[0]=r.x; a[1]=r.y; a[2]=r.z;
                a[3]=n.x; a[4]=n.y; a[5]=n.z;
            }
            b=(double)geomlib_vec3_dot(geomlib_vec3_sub(cp,q0),n);
            sum += b*b;
            for (k=0;k<6;++k) {
                atb[k]+=a[k]*b;
                for (l=0;l<6;++l) ata[k*6+l]+=a[k]*a[l];
            }
        }
        if (sum<=0.0) break;
        if (final_rms) *final_rms=(float)sqrt(sum/source->vertex_count);
        if (fabs(prev-(double)sqrt(sum/source->vertex_count))<p.tolerance) break;
        prev=sqrt(sum/source->vertex_count);
        if (geomlib_solve_gauss_d(6,ata,atb)!=GEOMLIB_SOLVER_SUCCESS) break;
        solve_rc=1;
        (void)solve_rc;
        {
            geomlib_vec3 w=geomlib_vec3_create((float)atb[0],(float)atb[1],(float)atb[2]);
            geomlib_vec3 t=geomlib_vec3_create((float)atb[3],(float)atb[4],(float)atb[5]);
            geomlib_mat4 U=rodrigues(w,t);
            T=geomlib_mat4_mul(U,T);
        }
    }
    geomlib_bvh_free(&bvh);
    *transform=T;
    return GEOMLIB_MESH_SUCCESS;
}
