/*
 * geomlib - Circle fitting (see geomlib_circle_fit.h)
 *
 * Clean-room algebraic least squares: solve for (a,b,c) in
 * x^2+y^2 + a x + b y + c = 0, then centre = (-a/2,-b/2).
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_circle_fit.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_fit_circle(const geomlib_vec2 *points, size_t count,
                       geomlib_circle_fit *out)
{
    double m[3][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    double x[3]={0,0,0};
    size_t i;
    int row,col,k;

    if (!out || !points || count<3) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<count;++i) {
        double px=points[i].x,py=points[i].y;
        double rhs=-(px*px+py*py);
        double r[4]={px,py,1.0,rhs};
        for (row=0;row<3;++row)
            for (col=0;col<4;++col)
                m[row][col]+=r[row]*r[col];
    }
    for (col=0;col<3;++col) {
        int piv=col;
        for (row=col+1;row<3;++row)
            if (fabs(m[row][col])>fabs(m[piv][col])) piv=row;
        if (fabs(m[piv][col])<1e-14) return GEOMLIB_MESH_ERROR_INVALID;
        if (piv!=col)
            for (k=0;k<4;++k) {
                double t=m[col][k]; m[col][k]=m[piv][k]; m[piv][k]=t;
            }
        for (row=col+1;row<3;++row) {
            double f=m[row][col]/m[col][col];
            for (k=0;k<4;++k) m[row][k]-=f*m[col][k];
        }
    }
    for (row=2;row>=0;--row) {
        double s=m[row][3];
        for (k=row+1;k<3;++k) s-=m[row][k]*x[k];
        x[row]=s/m[row][row];
    }
    {
        double cx=-0.5*x[0],cy=-0.5*x[1];
        double rr=cx*cx+cy*cy-x[2];
        double sum=0.0;
        if (!(rr>0.0)) return GEOMLIB_MESH_ERROR_INVALID;
        out->center=geomlib_vec2_create((float)cx,(float)cy);
        out->radius=(float)sqrt(rr);
        for (i=0;i<count;++i) {
            double d=hypot((double)points[i].x-cx,(double)points[i].y-cy)-sqrt(rr);
            sum+=d*d;
        }
        out->rms=(float)sqrt(sum/count);
    }
    return GEOMLIB_MESH_SUCCESS;
}
