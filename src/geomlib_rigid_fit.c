/*
 * geomlib - Rigid point-set fitting (see geomlib_rigid_fit.h)
 *
 * Clean-room Horn quaternion method with a 4x4 symmetric Jacobi eigensolve.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_rigid_fit.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Largest eigenvector of a symmetric 4x4 matrix by cyclic Jacobi rotations. */
static void largest_eigen4(double a[4][4], double out[4])
{
    double v[4][4];
    int sweep,p,q,i;
    for (i=0;i<4;++i)
        for (p=0;p<4;++p) v[i][p]=(i==p)?1.0:0.0;
    for (sweep=0;sweep<64;++sweep) {
        double off=0.0;
        for (p=0;p<4;++p) for (q=p+1;q<4;++q) off+=fabs(a[p][q]);
        if (off<1e-18) break;
        for (p=0;p<4;++p)
            for (q=p+1;q<4;++q) {
                double apq=a[p][q],tau,t,c,s;
                if (fabs(apq)<1e-22) continue;
                tau=(a[q][q]-a[p][p])/(2.0*apq);
                t=(tau>=0.0?1.0:-1.0)/(fabs(tau)+sqrt(1.0+tau*tau));
                c=1.0/sqrt(1.0+t*t); s=t*c;
                {
                    double app=a[p][p],aqq=a[q][q];
                    int k;
                    a[p][p]=c*c*app-2.0*s*c*apq+s*s*aqq;
                    a[q][q]=s*s*app+2.0*s*c*apq+c*c*aqq;
                    a[p][q]=a[q][p]=0.0;
                    for (k=0;k<4;++k) {
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
    {
        int best=0;
        for (i=1;i<4;++i) if (a[i][i]>a[best][best]) best=i;
        for (i=0;i<4;++i) out[i]=v[i][best];
    }
}

int geomlib_fit_rigid_points(const geomlib_vec3 *src, const geomlib_vec3 *dst,
                             size_t count, bool allow_scale,
                             geomlib_mat4 *transform, float *rms_out)
{
    double cs[3]={0,0,0},cd[3]={0,0,0};
    double sxx=0,sxy=0,sxz=0,syx=0,syy=0,syz=0,szx=0,szy=0,szz=0;
    double n[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    double q[4],scale=1.0;
    size_t i;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!transform || !src || !dst || count<3) return rc;
    for (i=0;i<count;++i) {
        cs[0]+=src[i].x; cs[1]+=src[i].y; cs[2]+=src[i].z;
        cd[0]+=dst[i].x; cd[1]+=dst[i].y; cd[2]+=dst[i].z;
    }
    for (i=0;i<3;++i) { cs[i]/=count; cd[i]/=count; }
    for (i=0;i<count;++i) {
        double px=(double)src[i].x-cs[0],py=(double)src[i].y-cs[1],pz=(double)src[i].z-cs[2];
        double qx=(double)dst[i].x-cd[0],qy=(double)dst[i].y-cd[1],qz=(double)dst[i].z-cd[2];
        sxx+=px*qx; sxy+=px*qy; sxz+=px*qz;
        syx+=py*qx; syy+=py*qy; syz+=py*qz;
        szx+=pz*qx; szy+=pz*qy; szz+=pz*qz;
    }
    n[0][0]=sxx+syy+szz;
    n[0][1]=syz-szy; n[0][2]=szx-sxz; n[0][3]=sxy-syx;
    n[1][1]=sxx-syy-szz; n[1][2]=sxy+syx; n[1][3]=szx+sxz;
    n[2][2]=-sxx+syy-szz; n[2][3]=syz+szy;
    n[3][3]=-sxx-syy+szz;
    for (i=0;i<4;++i) for (int j=i+1;j<4;++j) n[j][i]=n[i][j];
    largest_eigen4(n,q);
    {
        double w=q[0],x=q[1],y=q[2],z=q[3];
        double r[3][3];
        r[0][0]=1.0-2.0*(y*y+z*z); r[0][1]=2.0*(x*y-w*z); r[0][2]=2.0*(x*z+w*y);
        r[1][0]=2.0*(x*y+w*z); r[1][1]=1.0-2.0*(x*x+z*z); r[1][2]=2.0*(y*z-w*x);
        r[2][0]=2.0*(x*z-w*y); r[2][1]=2.0*(y*z+w*x); r[2][2]=1.0-2.0*(x*x+y*y);
        if (allow_scale) {
            double num=0.0,den=0.0;
            for (i=0;i<count;++i) {
                double px=(double)src[i].x-cs[0],py=(double)src[i].y-cs[1],pz=(double)src[i].z-cs[2];
                double qx=(double)dst[i].x-cd[0],qy=(double)dst[i].y-cd[1],qz=(double)dst[i].z-cd[2];
                double rx=r[0][0]*px+r[0][1]*py+r[0][2]*pz;
                double ry=r[1][0]*px+r[1][1]*py+r[1][2]*pz;
                double rz=r[2][0]*px+r[2][1]*py+r[2][2]*pz;
                num+=rx*qx+ry*qy+rz*qz;
                den+=px*px+py*py+pz*pz;
            }
            if (den>1e-18) scale=num/den;
        }
        for (int col=0;col<3;++col)
            for (int row=0;row<3;++row)
                transform->m[col*4+row]=(float)(scale*r[row][col]);
        transform->m[3]=0.0f;
        transform->m[7]=0.0f;
        transform->m[11]=0.0f;
        transform->m[15]=1.0f;
        transform->m[12]=(float)(cd[0]-scale*(r[0][0]*cs[0]+r[0][1]*cs[1]+r[0][2]*cs[2]));
        transform->m[13]=(float)(cd[1]-scale*(r[1][0]*cs[0]+r[1][1]*cs[1]+r[1][2]*cs[2]));
        transform->m[14]=(float)(cd[2]-scale*(r[2][0]*cs[0]+r[2][1]*cs[1]+r[2][2]*cs[2]));
    }
    if (rms_out) {
        double sum=0.0;
        for (i=0;i<count;++i) {
            geomlib_vec4 v=geomlib_mat4_mul_vec4(*transform,
                geomlib_vec4_create(src[i].x,src[i].y,src[i].z,1.0f));
            double dx=(double)v.x-dst[i].x,dy=(double)v.y-dst[i].y,dz=(double)v.z-dst[i].z;
            sum+=dx*dx+dy*dy+dz*dz;
        }
        *rms_out=(float)sqrt(sum/count);
    }
    return GEOMLIB_MESH_SUCCESS;
}
