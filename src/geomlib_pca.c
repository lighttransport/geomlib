/*
 * geomlib - PCA (see geomlib_pca.h)
 *
 * Clean-room covariance eigen decomposition by cyclic Jacobi rotations.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_pca.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_compute_pca(const geomlib_vec3 *points, size_t count,
                        geomlib_pca_result *out)
{
    double mean[3]={0,0,0};
    double cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    size_t i;
    int sweep,p,q,k,a,b,ord[3];

    if (!out||!points||count<2) return GEOMLIB_MESH_ERROR_INVALID;
    memset(out,0,sizeof(*out));
    for (i=0;i<count;++i) {
        mean[0]+=points[i].x; mean[1]+=points[i].y; mean[2]+=points[i].z;
    }
    for (k=0;k<3;++k) mean[k]/=count;
    for (i=0;i<count;++i) {
        double d[3]={(double)points[i].x-mean[0],
                     (double)points[i].y-mean[1],
                     (double)points[i].z-mean[2]};
        for (a=0;a<3;++a)
            for (b=0;b<3;++b) cov[a][b]+=d[a]*d[b];
    }
    for (a=0;a<3;++a) for (b=0;b<3;++b) cov[a][b]/=count;

    for (sweep=0;sweep<48;++sweep) {
        double off=0.0;
        for (p=0;p<2;++p) for (q=p+1;q<3;++q) off+=fabs(cov[p][q]);
        if (off<1e-20) break;
        for (p=0;p<2;++p)
            for (q=p+1;q<3;++q) {
                double apq=cov[p][q],tau,t,c,s;
                if (fabs(apq)<1e-22) continue;
                tau=(cov[q][q]-cov[p][p])/(2.0*apq);
                t=(tau>=0.0?1.0:-1.0)/(fabs(tau)+sqrt(1.0+tau*tau));
                c=1.0/sqrt(1.0+t*t); s=t*c;
                {
                    double app=cov[p][p],aqq=cov[q][q],apq2=apq;
                    int r;
                    cov[p][p]=c*c*app-2.0*s*c*apq2+s*s*aqq;
                    cov[q][q]=s*s*app+2.0*s*c*apq2+c*c*aqq;
                    cov[p][q]=cov[q][p]=0.0;
                    for (r=0;r<3;++r) {
                        if (r!=p&&r!=q) {
                            double akp=cov[r][p],akq=cov[r][q];
                            cov[r][p]=cov[p][r]=c*akp-s*akq;
                            cov[r][q]=cov[q][r]=s*akp+c*akq;
                        }
                        {
                            double vkp=v[r][p],vkq=v[r][q];
                            v[r][p]=c*vkp-s*vkq;
                            v[r][q]=s*vkp+c*vkq;
                        }
                    }
                }
            }
    }
    ord[0]=0; ord[1]=1; ord[2]=2;
    for (p=0;p<3;++p)
        for (q=p+1;q<3;++q)
            if (cov[ord[q]][ord[q]]>cov[ord[p]][ord[p]]) {
                int t=ord[p]; ord[p]=ord[q]; ord[q]=t;
            }
    for (k=0;k<3;++k) {
        int a=ord[k];
        out->eigenvalues[k]=cov[a][a];
        out->eigenvectors[k]=geomlib_vec3_create((float)v[0][a],(float)v[1][a],(float)v[2][a]);
    }
    out->mean=geomlib_vec3_create((float)mean[0],(float)mean[1],(float)mean[2]);
    return GEOMLIB_MESH_SUCCESS;
}
