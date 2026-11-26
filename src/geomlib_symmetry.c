/*
 * geomlib - Symmetry detection (see geomlib_symmetry.h)
 *
 * Clean-room principal-axis plane testing with BVH reflection matching.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_symmetry.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Eigenvectors of symmetric 3x3 (columns of out). */
static void eigen3(double a[3][3], double out[3][3])
{
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    int sweep,p,q;
    for (sweep=0;sweep<32;++sweep) {
        double off=fabs(a[0][1])+fabs(a[0][2])+fabs(a[1][2]);
        if (off<1e-20) break;
        for (p=0;p<2;++p)
            for (q=p+1;q<3;++q) {
                double apq=a[p][q],tau,t,c,s;
                if (fabs(apq)<1e-22) continue;
                tau=(a[q][q]-a[p][p])/(2.0*apq);
                t=(tau>=0.0?1.0:-1.0)/(fabs(tau)+sqrt(1.0+tau*tau));
                c=1.0/sqrt(1.0+t*t); s=t*c;
                {
                    double app=a[p][p],aqq=a[q][q],apq2=apq;
                    int k;
                    a[p][p]=c*c*app-2.0*s*c*apq2+s*s*aqq;
                    a[q][q]=s*s*app+2.0*s*c*apq2+c*c*aqq;
                    a[p][q]=a[q][p]=0.0;
                    for (k=0;k<3;++k) {
                        if (k!=p && k!=q) {
                            double akp=a[k][p],akq=a[k][q];
                            a[k][p]=a[p][k]=c*akp-s*akq;
                            a[k][q]=a[q][k]=s*akp+c*akq;
                        }
                        {
                            double vkp=v[k][p],vkq=v[k][q];
                            v[k][p]=c*vkp-s*vkq;
                            v[k][q]=s*vkp+c*vkq;
                        }
                    }
                }
            }
    }
    for (p=0;p<3;++p)
        for (q=0;q<3;++q) out[p][q]=v[p][q];
}

int geomlib_detect_symmetry(const geomlib_trimesh *mesh,
                            geomlib_symmetry_result *out)
{
    geomlib_bvh bvh;
    double c[3]={0,0,0},cov[3][3]={{0,0,0},{0,0,0},{0,0,0}},ev[3][3];
    size_t i,k;
    int axis,best_axis=0;
    double best=DBL_MAX;
    double best_d[3]={0,0,1};

    if (!out || !mesh || !mesh->vertices || !mesh->triangle_count ||
        !mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_bvh_build(&bvh,mesh,GEOMLIB_BVH_QUALITY_MEDIUM)!=
        GEOMLIB_MESH_SUCCESS) return GEOMLIB_MESH_ERROR_INVALID;

    for (i=0;i<mesh->vertex_count;++i) {
        c[0]+=mesh->vertices[i].x;
        c[1]+=mesh->vertices[i].y;
        c[2]+=mesh->vertices[i].z;
    }
    for (k=0;k<3;++k) c[k]/=mesh->vertex_count;
    for (i=0;i<mesh->vertex_count;++i) {
        double d[3]={(double)mesh->vertices[i].x-c[0],
                     (double)mesh->vertices[i].y-c[1],
                     (double)mesh->vertices[i].z-c[2]};
        for (k=0;k<3;++k) {
            double v[3]={d[0],d[1],d[2]};
            int j;
            for (j=0;j<3;++j) cov[k][j]+=v[k]*v[j];
        }
    }
    for (k=0;k<3;++k) for (i=0;i<3;++i) cov[k][i]/=mesh->vertex_count;
    eigen3(cov,ev);

    for (axis=0;axis<3;++axis) {
        double n[3]={ev[0][axis],ev[1][axis],ev[2][axis]};
        double len=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        double sum=0.0;
        if (!(len>1e-12)) continue;
        for (k=0;k<3;++k) n[k]/=len;
        for (i=0;i<mesh->vertex_count;++i) {
            double d=(double)mesh->vertices[i].x-c[0];
            double e=(double)mesh->vertices[i].y-c[1];
            double f=(double)mesh->vertices[i].z-c[2];
            double t=2.0*(d*n[0]+e*n[1]+f*n[2]);
            geomlib_vec3 q=geomlib_vec3_create(
                (float)(mesh->vertices[i].x-t*n[0]),
                (float)(mesh->vertices[i].y-t*n[1]),
                (float)(mesh->vertices[i].z-t*n[2]));
            geomlib_vec3 cp; float dist=FLT_MAX;
            if (geomlib_bvh_closest_point(&bvh,q,&cp,&dist,NULL) && isfinite(dist))
                sum+=(double)dist*dist;
        }
        {
            double rms=sqrt(sum/mesh->vertex_count);
            if (rms<best) { best=rms; best_axis=axis;
                            for (k=0;k<3;++k) best_d[k]=n[k]; }
        }
    }
    geomlib_bvh_free(&bvh);
    out->plane_point=geomlib_vec3_create((float)c[0],(float)c[1],(float)c[2]);
    out->plane_normal=geomlib_vec3_create((float)best_d[0],(float)best_d[1],
                                          (float)best_d[2]);
    out->rms_error=(float)best;
    (void)best_axis;
    return GEOMLIB_MESH_SUCCESS;
}
