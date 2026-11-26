/*
 * geomlib - Point-cloud curvature (see geomlib_point_curvature.h)
 *
 * Clean-room local PCA with cyclic Jacobi eigen decomposition.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_point_curvature.h"
#include "geomlib/geomlib_kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Smallest eigenvalue and eigenvector of a symmetric 3x3 matrix. */
static double smallest_eigen(double a[3][3], double out[3])
{
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    int sweep,p,q,m;
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
    m=0;
    if (a[1][1]<a[m][m]) m=1;
    if (a[2][2]<a[m][m]) m=2;
    out[0]=v[0][m]; out[1]=v[1][m]; out[2]=v[2][m];
    return a[m][m];
}

int geomlib_pointcloud_curvature(const geomlib_pointcloud *cloud, int k,
                                 float *curvature, geomlib_vec3 *normals)
{
    geomlib_kdtree *tree=NULL;
    uint32_t *idx=NULL;
    size_t n, found, i, j;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!cloud || !cloud->points || !cloud->point_count ||
        (!curvature && !normals)) return rc;
    if (k<=0) k=16;
    n=(size_t)k;
    if (n>cloud->point_count) n=cloud->point_count;
    if (n<3) return rc;
    tree=geomlib_kdtree_build(cloud->points,cloud->point_count);
    idx=(uint32_t*)malloc(n*sizeof(uint32_t));
    if (!tree||!idx) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }

    for (i=0;i<cloud->point_count;++i) {
        double cx=0.0,cy=0.0,cz=0.0;
        double cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        double ev[3],sum;
        found=geomlib_kdtree_knn(tree,cloud->points[i],n,idx,NULL);
        if (found<3) continue;
        for (j=0;j<found;++j) {
            geomlib_vec3 p=cloud->points[idx[j]];
            cx+=p.x; cy+=p.y; cz+=p.z;
        }
        cx/=found; cy/=found; cz/=found;
        for (j=0;j<found;++j) {
            geomlib_vec3 p=cloud->points[idx[j]];
            double dx=(double)p.x-cx,dy=(double)p.y-cy,dz=(double)p.z-cz;
            cov[0][0]+=dx*dx; cov[0][1]+=dx*dy; cov[0][2]+=dx*dz;
            cov[1][1]+=dy*dy; cov[1][2]+=dy*dz; cov[2][2]+=dz*dz;
        }
        cov[1][0]=cov[0][1]; cov[2][0]=cov[0][2]; cov[2][1]=cov[1][2];
        smallest_eigen(cov,ev);
        sum=cov[0][0]+cov[1][1]+cov[2][2];
        if (curvature)
            curvature[i]=(float)(sum>1e-18 ? ev[0]/sum : 0.0);
        if (normals) {
            geomlib_vec3 v=geomlib_vec3_create((float)ev[0],(float)ev[1],(float)ev[2]);
            float len=geomlib_vec3_length(v);
            normals[i]=len>1e-12f ? geomlib_vec3_mul(v,1.0f/len)
                                  : geomlib_vec3_create(0,0,1);
        }
    }
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_kdtree_free(tree);
    free(idx);
    return rc;
}
