/*
 * geomlib - Plane fitting (see geomlib_plane_fit.h)
 *
 * Clean-room PCA with a cyclic Jacobi 3x3 eigensolve.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_plane_fit.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void smallest_eigen(double a[3][3], double out[3])
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
}

int geomlib_fit_plane(const geomlib_vec3 *points, size_t count,
                      geomlib_plane_fit *out)
{
    double c[3]={0,0,0};
    double cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double n[3];
    size_t i;
    double sum=0.0;

    if (!out || !points || count<3) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<count;++i) {
        c[0]+=points[i].x; c[1]+=points[i].y; c[2]+=points[i].z;
    }
    for (i=0;i<3;++i) c[i]/=count;
    for (i=0;i<count;++i) {
        double d[3]={(double)points[i].x-c[0],
                     (double)points[i].y-c[1],
                     (double)points[i].z-c[2]};
        int j,k;
        for (j=0;j<3;++j)
            for (k=0;k<3;++k) cov[j][k]+=d[j]*d[k];
    }
    for (i=0;i<9;++i) ((double*)cov)[i]/=count;
    smallest_eigen(cov,n);
    {
        double len=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (!(len>1e-12)) return GEOMLIB_MESH_ERROR_INVALID;
        n[0]/=len; n[1]/=len; n[2]/=len;
    }
    for (i=0;i<count;++i) {
        double d=((double)points[i].x-c[0])*n[0]
               + ((double)points[i].y-c[1])*n[1]
               + ((double)points[i].z-c[2])*n[2];
        sum+=d*d;
    }
    out->point=geomlib_vec3_create((float)c[0],(float)c[1],(float)c[2]);
    out->normal=geomlib_vec3_create((float)n[0],(float)n[1],(float)n[2]);
    out->rms=(float)sqrt(sum/count);
    return GEOMLIB_MESH_SUCCESS;
}
