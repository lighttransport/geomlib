/*
 * geomlib - Sphere fitting (see geomlib_sphere_fit.h)
 *
 * Clean-room algebraic least squares: solve for (a,b,c,d) in
 * x^2+y^2+z^2 + a x + b y + c z + d = 0, centre = (-a/2,-b/2,-c/2).
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sphere_fit.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_fit_sphere(const geomlib_vec3 *points, size_t count,
                       geomlib_sphere_fit *out)
{
    double m[4][5]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    double x[4]={0,0,0,0};
    size_t i;
    int row,col,k;

    if (!out || !points || count<4) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<count;++i) {
        double px=points[i].x,py=points[i].y,pz=points[i].z;
        double rhs=-(px*px+py*py+pz*pz);
        double r[5]={px,py,pz,1.0,rhs};
        for (row=0;row<4;++row)
            for (col=0;col<5;++col)
                m[row][col]+=r[row]*r[col];
    }
    for (col=0;col<4;++col) {
        int piv=col;
        for (row=col+1;row<4;++row)
            if (fabs(m[row][col])>fabs(m[piv][col])) piv=row;
        if (fabs(m[piv][col])<1e-14) return GEOMLIB_MESH_ERROR_INVALID;
        if (piv!=col)
            for (k=0;k<5;++k) {
                double t=m[col][k]; m[col][k]=m[piv][k]; m[piv][k]=t;
            }
        for (row=col+1;row<4;++row) {
            double f=m[row][col]/m[col][col];
            for (k=0;k<5;++k) m[row][k]-=f*m[col][k];
        }
    }
    for (row=3;row>=0;--row) {
        double s=m[row][4];
        for (k=row+1;k<4;++k) s-=m[row][k]*x[k];
        x[row]=s/m[row][row];
    }
    {
        double cx=-0.5*x[0],cy=-0.5*x[1],cz=-0.5*x[2];
        double rr=cx*cx+cy*cy+cz*cz-x[3];
        double sum=0.0;
        if (!(rr>0.0)) return GEOMLIB_MESH_ERROR_INVALID;
        out->center=geomlib_vec3_create((float)cx,(float)cy,(float)cz);
        out->radius=(float)sqrt(rr);
        for (i=0;i<count;++i) {
            double d=geomlib_vec3_length(geomlib_vec3_sub(points[i],out->center))-sqrt(rr);
            sum+=d*d;
        }
        out->rms=(float)sqrt(sum/count);
    }
    return GEOMLIB_MESH_SUCCESS;
}
