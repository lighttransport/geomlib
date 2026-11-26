/*
 * geomlib - Minimum enclosing sphere (see geomlib_miniball.h)
 *
 * Clean-room Welzl randomized incremental algorithm with exact boundary
 * sphere construction for up to four support points.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_miniball.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { geomlib_vec3 c; float r; int valid; } ball;

static float dist(const geomlib_vec3 *a, const geomlib_vec3 *b)
{
    return geomlib_vec3_length(geomlib_vec3_sub(*a,*b));
}

static int contains(const ball *b, geomlib_vec3 p, float eps)
{
    if (!b->valid) return 0;
    return geomlib_vec3_length(geomlib_vec3_sub(p,b->c)) <= b->r+eps;
}

/* Sphere through up to 4 points; returns invalid when degenerate. */
static ball sphere_from(const geomlib_vec3 *p, size_t n)
{
    ball b;
    b.valid=0; b.r=-1.0f; b.c=geomlib_vec3_create(0,0,0);
    if (n==1) { b.c=p[0]; b.r=0.0f; b.valid=1; return b; }
    if (n==2) {
        b.c=geomlib_vec3_mul(geomlib_vec3_add(p[0],p[1]),0.5f);
        b.r=0.5f*dist(&p[0],&p[1]); b.valid=b.r>=0.0f; return b;
    }
    {
        /* Solve for centre in the affine span of the support points. */
        geomlib_vec3 o=p[0];
        double a[3][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0}};
        double x[3]={0,0,0};
        size_t i,j,k;
        int m=(int)n-1;
        for (i=0;i<(size_t)m;++i) {
            geomlib_vec3 d=geomlib_vec3_sub(p[i+1],o);
            a[i][0]=d.x; a[i][1]=d.y; a[i][2]=d.z;
            a[i][3]=0.5*((double)d.x*d.x+(double)d.y*d.y+(double)d.z*d.z);
        }
        for (i=0;i<(size_t)m;++i) {
            int piv=(int)i;
            for (j=(size_t)i+1;j<(size_t)m;++j)
                if (fabs(a[j][i])>fabs(a[piv][i])) piv=(int)j;
            if (fabs(a[piv][i])<1e-12) return b;
            if (piv!=(int)i)
                for (k=0;k<4;++k) {
                    double t=a[i][k]; a[i][k]=a[piv][k]; a[piv][k]=t;
                }
            for (j=(size_t)i+1;j<(size_t)m;++j) {
                double f=a[j][i]/a[i][i];
                for (k=0;k<4;++k) a[j][k]-=f*a[i][k];
            }
        }
        for (i=(size_t)m;i-- > 0;) {
            double s=a[i][3];
            for (j=i+1;j<(size_t)m;++j) s-=a[i][j]*x[j];
            x[i]=s/a[i][i];
        }
        {
            geomlib_vec3 c=o;
            geomlib_vec3 v=geomlib_vec3_create((float)x[0],(float)x[1],(float)x[2]);
            c=geomlib_vec3_add(c,v);
            b.c=c;
            b.r=dist(&c,&p[0]);
            b.valid=1;
        }
    }
    return b;
}

/* Smallest ball enclosing the boundary set: try every supporting subset. */
static ball trivial(const geomlib_vec3 *p, size_t n)
{
    ball best;
    unsigned mask;
    size_t i;
    best.valid=0; best.r=FLT_MAX; best.c=geomlib_vec3_create(0,0,0);
    for (mask=1;mask<(1u<<n);++mask) {
        geomlib_vec3 sub[4];
        size_t m=0,cnt=0;
        ball b;
        for (i=0;i<n;++i) if (mask&(1u<<i)) sub[m++]=p[i];
        b=sphere_from(sub,m);
        if (!b.valid) continue;
        for (i=0;i<n;++i)
            if (!contains(&b,p[i],1e-5f)) { cnt=1; break; }
        if (!cnt && b.r<best.r) best=b;
    }
    return best;
}

static ball welzl(const geomlib_vec3 *p, size_t n,
                  const geomlib_vec3 *r, size_t nr)
{
    ball d;
    size_t i;
    if (n==0 || nr==4) return trivial(r,nr);
    d=welzl(p,n-1,r,nr);
    if (contains(&d,p[n-1],1e-5f)) return d;
    {
        geomlib_vec3 nr2[4];
        for (i=0;i<nr;++i) nr2[i]=r[i];
        nr2[nr]=p[n-1];
        return welzl(p,n-1,nr2,nr+1);
    }
}

int geomlib_minimum_enclosing_sphere(const geomlib_vec3 *points, size_t count,
                                     geomlib_bounding_sphere *out)
{
    ball b;
    if (!out || !points || !count) return GEOMLIB_MESH_ERROR_INVALID;
    if (count>512) return GEOMLIB_MESH_ERROR_INVALID;
    b=welzl(points,count,NULL,0);
    if (!b.valid) return GEOMLIB_MESH_ERROR_INVALID;
    out->center=b.c; out->radius=b.r;
    return GEOMLIB_MESH_SUCCESS;
}
