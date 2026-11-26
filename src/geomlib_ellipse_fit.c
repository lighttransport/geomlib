/*
 * geomlib - Ellipse fitting (see geomlib_ellipse_fit.h)
 *
 * Clean-room nonlinear least squares on (centre, semi-axes, angle) with PCA
 * initialization and backtracking gradient descent.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_ellipse_fit.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double cx,cy,A,B,t; } ellipse_params;

static double energy(const geomlib_vec2 *p, size_t n, ellipse_params e)
{
    double sum=0.0;
    size_t i;
    double cs=cos(e.t),sn=sin(e.t);
    for (i=0;i<n;++i) {
        double dx=(double)p[i].x-e.cx,dy=(double)p[i].y-e.cy;
        double u=dx*cs+dy*sn, v=-dx*sn+dy*cs;
        double r=sqrt((u*u)/(e.A*e.A)+(v*v)/(e.B*e.B))-1.0;
        sum+=r*r;
    }
    return sum;
}

static void gradient(const geomlib_vec2 *p, size_t n, ellipse_params e,
                     double g[5])
{
    const double h=1e-5;
    int k;
    for (k=0;k<5;++k) {
        ellipse_params a=e,b=e;
        switch (k) {
            case 0: a.cx+=h; b.cx-=h; break;
            case 1: a.cy+=h; b.cy-=h; break;
            case 2: a.A+=h;  b.A-=h;  break;
            case 3: a.B+=h;  b.B-=h;  break;
            default: a.t+=h; b.t-=h; break;
        }
        if (a.A<=1e-6) a.A=1e-6;
        if (b.A<=1e-6) b.A=1e-6;
        if (a.B<=1e-6) a.B=1e-6;
        if (b.B<=1e-6) b.B=1e-6;
        g[k]=(energy(p,n,a)-energy(p,n,b))/(2.0*h);
    }
}

int geomlib_fit_ellipse(const geomlib_vec2 *points, size_t count,
                        geomlib_ellipse_fit *out)
{
    ellipse_params e;
    double cx=0.0,cy=0.0,cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    size_t i;
    int sweep,p,q,iter;

    if (!out || !points || count<5) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<count;++i) { cx+=points[i].x; cy+=points[i].y; }
    cx/=count; cy/=count;
    for (i=0;i<count;++i) {
        double dx=(double)points[i].x-cx,dy=(double)points[i].y-cy;
        cov[0][0]+=dx*dx; cov[0][1]+=dx*dy; cov[1][1]+=dy*dy;
    }
    cov[0][0]/=count; cov[0][1]/=count; cov[1][1]/=count;
    cov[1][0]=cov[0][1];
    /* 2x2 Jacobi rotation for the covariance eigenvectors. */
    {
        double th=0.5*atan2(2.0*cov[0][1],cov[0][0]-cov[1][1]);
        double cs=cos(th),sn=sin(th);
        double l1=cs*cs*cov[0][0]+2.0*cs*sn*cov[0][1]+sn*sn*cov[1][1];
        double l2=sn*sn*cov[0][0]-2.0*cs*sn*cov[0][1]+cs*cs*cov[1][1];
        e.cx=cx; e.cy=cy;
        e.A=2.0*sqrt(fmax(l1,1e-9));
        e.B=2.0*sqrt(fmax(l2,1e-9));
        e.t=th;
        if (e.A<e.B) { double t=e.A; e.A=e.B; e.B=t; e.t+=M_PI*0.5; }
    }
    (void)sweep; (void)p; (void)q; (void)v;

    for (iter=0;iter<200;++iter) {
        double g[5],step=0.1;
        double e0=energy(points,count,e);
        gradient(points,count,e,g);
        {
            double mx=0.0;
            for (p=0;p<5;++p) if (fabs(g[p])>mx) mx=fabs(g[p]);
            if (mx<1e-12) break;
            step=0.1/mx;
        }
        for (q=0;q<20;++q) {
            ellipse_params t=e;
            t.cx-=step*g[0]; t.cy-=step*g[1];
            t.A-=step*g[2];  t.B-=step*g[3];  t.t-=step*g[4];
            if (t.A>1e-6 && t.B>1e-6 && energy(points,count,t)<e0) {
                e=t; break;
            }
            step*=0.5;
        }
        if (q==20) break;
    }
    if (e.B>e.A) { double t=e.A; e.A=e.B; e.B=t; e.t+=M_PI*0.5; }
    out->center=geomlib_vec2_create((float)e.cx,(float)e.cy);
    out->axis_a=(float)e.A;
    out->axis_b=(float)e.B;
    out->angle=(float)fmod(e.t,M_PI);
    {
        double cs=cos(e.t),sn=sin(e.t),sum=0.0;
        for (i=0;i<count;++i) {
            double dx=(double)points[i].x-e.cx,dy=(double)points[i].y-e.cy;
            double u=dx*cs+dy*sn,v=-dx*sn+dy*cs;
            double rr=sqrt((u*u)/(e.A*e.A)+(v*v)/(e.B*e.B))-1.0;
            sum+=rr*rr;
        }
        out->rms=(float)sqrt(sum/count);
    }
    return GEOMLIB_MESH_SUCCESS;
}
