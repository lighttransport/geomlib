/*
 * geomlib - Cylinder fitting (see geomlib_cylinder_fit.h)
 *
 * Clean-room nonlinear least squares on (axis, point, radius) with gradient
 * descent and backtracking. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_cylinder_fit.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { geomlib_vec3 a,p; float r; } cyl_params;

static double energy(const geomlib_vec3 *pts, size_t n, cyl_params c)
{
    double sum=0.0;
    size_t i;
    for (i=0;i<n;++i) {
        geomlib_vec3 d=geomlib_vec3_sub(pts[i],c.p);
        float t=geomlib_vec3_dot(d,c.a);
        geomlib_vec3 perp=geomlib_vec3_sub(d,geomlib_vec3_mul(c.a,t));
        double r=geomlib_vec3_length(perp)-(double)c.r;
        sum+=r*r;
    }
    return sum;
}

static void normalize_axis(geomlib_vec3 *a)
{
    float l=geomlib_vec3_length(*a);
    if (l>1e-12f) *a=geomlib_vec3_mul(*a,1.0f/l);
    else *a=geomlib_vec3_create(0,0,1);
}

int geomlib_fit_cylinder(const geomlib_vec3 *points, size_t count,
                         geomlib_vec3 initial_axis,
                         geomlib_cylinder_fit *out)
{
    cyl_params c;
    double cx=0,cy=0,cz=0;
    size_t i;
    int iter,k,q;

    if (!out || !points || count<5) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<count;++i) { cx+=points[i].x; cy+=points[i].y; cz+=points[i].z; }
    c.p=geomlib_vec3_create((float)(cx/count),(float)(cy/count),(float)(cz/count));
    c.a=initial_axis;
    normalize_axis(&c.a);
    {
        double sum=0.0;
        for (i=0;i<count;++i) {
            geomlib_vec3 d=geomlib_vec3_sub(points[i],c.p);
            float t=geomlib_vec3_dot(d,c.a);
            geomlib_vec3 perp=geomlib_vec3_sub(d,geomlib_vec3_mul(c.a,t));
            sum+=geomlib_vec3_length(perp);
        }
        c.r=(float)(sum/count);
        if (!(c.r>1e-9f)) return GEOMLIB_MESH_ERROR_INVALID;
    }

    for (iter=0;iter<200;++iter) {
        double g[7],step;
        double e0=energy(points,count,c);
        const float h=1e-4f;
        for (k=0;k<7;++k) {
            cyl_params a=c,b=c;
            switch (k) {
                case 0: a.a.x+=h; b.a.x-=h; break;
                case 1: a.a.y+=h; b.a.y-=h; break;
                case 2: a.a.z+=h; b.a.z-=h; break;
                case 3: a.p.x+=h; b.p.x-=h; break;
                case 4: a.p.y+=h; b.p.y-=h; break;
                case 5: a.p.z+=h; b.p.z-=h; break;
                default: a.r+=h; b.r-=h; break;
            }
            normalize_axis(&a.a); normalize_axis(&b.a);
            if (a.r<=1e-6f) a.r=1e-6f;
            if (b.r<=1e-6f) b.r=1e-6f;
            g[k]=(energy(points,count,a)-energy(points,count,b))/(2.0*h);
        }
        {
            double mx=0.0;
            for (k=0;k<7;++k) if (fabs(g[k])>mx) mx=fabs(g[k]);
            if (mx<1e-12) break;
            step=0.05/mx;
        }
        for (q=0;q<20;++q) {
            cyl_params t=c;
            t.a.x-=step*g[0]; t.a.y-=step*g[1]; t.a.z-=step*g[2];
            t.p.x-=step*g[3]; t.p.y-=step*g[4]; t.p.z-=step*g[5];
            t.r-=(float)(step*g[6]);
            normalize_axis(&t.a);
            if (t.r>1e-6f && energy(points,count,t)<e0) { c=t; break; }
            step*=0.5;
        }
        if (q==20) break;
    }
    normalize_axis(&c.a);
    out->axis=c.a;
    out->point=c.p;
    out->radius=c.r;
    out->rms=(float)sqrt(energy(points,count,c)/count);
    return GEOMLIB_MESH_SUCCESS;
}
