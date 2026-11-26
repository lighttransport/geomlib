/*
 * geomlib - 2D Delaunay refinement (see geomlib_refine2d.h)
 *
 * Clean-room Chew/Ruppert style refinement: insert circumcentres of bad
 * triangles when admissible, otherwise split their longest edge (including
 * constrained boundary edges).
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_refine2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

geomlib_refine2d_params geomlib_refine2d_params_default(void)
{
    geomlib_refine2d_params p;
    p.min_angle_deg=25.0f;
    p.max_points=2000;
    return p;
}

static double orient(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x-a.x)*((double)c.y-a.y)
         - ((double)b.y-a.y)*((double)c.x-a.x);
}

static double corner_angle(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    double ux=(double)b.x-a.x, uy=(double)b.y-a.y;
    double vx=(double)c.x-a.x, vy=(double)c.y-a.y;
    double nu=hypot(ux,uy), nv=hypot(vx,vy), cs;
    if (!(nu>1e-14 && nv>1e-14)) return 0.0;
    cs=(ux*vx+uy*vy)/(nu*nv);
    if (cs>1.0) cs=1.0;
    if (cs<-1.0) cs=-1.0;
    return acos(cs);
}

float geomlib_cdt2d_min_angle(const geomlib_cdt2d *c)
{
    size_t i;
    double mn=DBL_MAX;
    if (!c || !c->triangles || !c->triangle_count) return 0.0f;
    for (i=0;i<c->triangle_count;++i) {
        geomlib_vec2 a=c->vertices[c->triangles[i*3]];
        geomlib_vec2 b=c->vertices[c->triangles[i*3+1]];
        geomlib_vec2 d=c->vertices[c->triangles[i*3+2]];
        double x=corner_angle(a,b,d), y=corner_angle(b,d,a), z=corner_angle(d,a,b);
        double t=x<y?(x<z?x:z):(y<z?y:z);
        if (t<mn) mn=t;
    }
    return (float)(mn*180.0/M_PI);
}

static int in_polygon(geomlib_vec2 q, const geomlib_vec2 *p, size_t n, double eps)
{
    size_t i;
    int inside=0;
    for (i=0;i<n;++i) {
        size_t j=(i+1)%n;
        if (fabs(orient(p[i],p[j],q))<eps) return 0;
        if (((p[i].y>q.y)!=(p[j].y>q.y)) &&
            ((double)q.x < ((double)p[j].x-p[i].x)*((double)q.y-p[i].y)/
                           ((double)p[j].y-p[i].y)+p[i].x))
            inside=!inside;
    }
    return inside;
}

static geomlib_vec2 circumcenter(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    double ax=a.x,ay=a.y,bx=b.x,by=b.y,cx=c.x,cy=c.y;
    double d=2.0*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by));
    if (fabs(d)<1e-18) return geomlib_vec2_create((a.x+b.x+c.x)/3.0f,
                                                   (a.y+b.y+c.y)/3.0f);
    {
        double a2=ax*ax+ay*ay, b2=bx*bx+by*by, c2=cx*cx+cy*cy;
        double x=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d;
        double y=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d;
        return geomlib_vec2_create((float)x,(float)y);
    }
}

/* Split a constraint (a,b) at new vertex m: replace it with (a,m),(m,b). */
static void split_constraint(uint32_t *cons, size_t *nc,
                             uint32_t a, uint32_t b, uint32_t m)
{
    size_t i;
    for (i=0;i+1<*nc;i+=2) {
        if ((cons[i]==a&&cons[i+1]==b)||(cons[i]==b&&cons[i+1]==a)) {
            cons[i+1]=m;
            cons[*nc]=m; cons[*nc+1]=b;
            (*nc)+=2;
            return;
        }
    }
}

int geomlib_refine_polygon2d(const geomlib_vec2 *polygon, size_t n,
                             const geomlib_refine2d_params *params,
                             geomlib_cdt2d *out)
{
    geomlib_refine2d_params p=params ? *params : geomlib_refine2d_params_default();
    geomlib_vec2 *pts=NULL;
    uint32_t *cons=NULL;
    size_t np=0, nc=0, i, iter;
    double target, extent=0.0;
    geomlib_cdt2d cdt;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !polygon || n<3) return rc;
    if (!(p.min_angle_deg>0.0f && p.min_angle_deg<60.0f)) p.min_angle_deg=25.0f;
    if (p.max_points<1) p.max_points=1;
    if (p.max_points > (1<<22)) p.max_points=1<<22;
    target=p.min_angle_deg*M_PI/180.0;
    for (i=0;i<n;++i) extent=fmax(extent,fmax(fabs(polygon[i].x),fabs(polygon[i].y)));
    if (!(extent>1e-9)) return rc;

    pts=(geomlib_vec2*)malloc((n+p.max_points)*sizeof(geomlib_vec2));
    cons=(uint32_t*)malloc((2*n+4*p.max_points)*sizeof(uint32_t));
    if (!pts||!cons) goto cleanup;
    for (i=0;i<n;++i) pts[np++]=polygon[i];
    for (i=0;i<n;++i) {
        cons[nc++]=(uint32_t)i;
        cons[nc++]=(uint32_t)((i+1)%n);
    }

    for (iter=0;iter<(size_t)p.max_points;++iter) {
        size_t t, worst=(size_t)-1;
        double wa=target, rad=extent*1e-8;
        cdt=geomlib_cdt2d_create();
        if (geomlib_cdt2d_build(pts,np,cons,nc/2,&cdt)!=GEOMLIB_MESH_SUCCESS)
            goto cleanup;
        for (t=0;t<cdt.triangle_count;++t) {
            geomlib_vec2 a=pts[cdt.triangles[t*3]],b=pts[cdt.triangles[t*3+1]],
                         c=pts[cdt.triangles[t*3+2]];
            double x=corner_angle(a,b,c),y=corner_angle(b,c,a),z=corner_angle(c,a,b);
            double mn=x<y?(x<z?x:z):(y<z?y:z);
            if (mn<wa) { wa=mn; worst=t; }
        }
        if (worst==(size_t)-1) { *out=cdt; goto done; }
        {
            uint32_t a=cdt.triangles[worst*3],b=cdt.triangles[worst*3+1],
                     c=cdt.triangles[worst*3+2];
            geomlib_vec2 cc=circumcenter(pts[a],pts[b],pts[c]);
            int too_close=0;
            for (i=0;i<np;++i)
                if (hypot((double)cc.x-pts[i].x,(double)cc.y-pts[i].y)<rad) {
                    too_close=1; break;
                }
            if (!too_close && in_polygon(cc,polygon,n,rad)) {
                pts[np++]=cc;
                geomlib_cdt2d_free(&cdt);
                continue;
            }
            /* Circumcentre inadmissible: split the longest edge of the
             * offending triangle at its midpoint. */
            {
                uint32_t u=a,v=b;
                double best=-1.0;
                int k;
                uint32_t vs[3]={a,b,c};
                for (k=0;k<3;++k) {
                    uint32_t x=vs[k],y=vs[(k+1)%3];
                    double len=hypot((double)pts[x].x-pts[y].x,
                                     (double)pts[x].y-pts[y].y);
                    if (len>best) { best=len; u=x; v=y; }
                }
                if (np<(size_t)p.max_points) {
                    uint32_t m=(uint32_t)np;
                    pts[np++]=geomlib_vec2_create((pts[u].x+pts[v].x)*0.5f,
                                                  (pts[u].y+pts[v].y)*0.5f);
                    split_constraint(cons,&nc,u,v,m);
                    geomlib_cdt2d_free(&cdt);
                    continue;
                }
            }
        }
        *out=cdt;
        goto done;
    }
    cdt=geomlib_cdt2d_create();
    geomlib_cdt2d_build(pts,np,cons,nc/2,&cdt);

done:
    rc=GEOMLIB_MESH_SUCCESS;
    *out=cdt;
cleanup:
    free(pts); free(cons);
    return rc;
}
