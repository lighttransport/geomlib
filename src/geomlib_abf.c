/*
 * geomlib - Angle-based flattening (see geomlib_abf.h)
 *
 * Clean-room nonlinear angle-distortion minimization. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_abf.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_abf_params geomlib_abf_params_default(void)
{
    geomlib_abf_params p;
    p.iterations=60;
    p.step=0.05f;
    p.tolerance=1e-3f;
    return p;
}

static double angle3(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c)
{
    geomlib_vec3 u=geomlib_vec3_sub(b,a), v=geomlib_vec3_sub(c,a);
    double nu=geomlib_vec3_length(u), nv=geomlib_vec3_length(v);
    if (!(nu>1e-12 && nv>1e-12)) return 0.0;
    {
        double cs=geomlib_vec3_dot(u,v)/(nu*nv);
        if (cs>1.0) cs=1.0;
        if (cs<-1.0) cs=-1.0;
        return acos(cs);
    }
}

static double angle2(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    double ux=(double)b.x-a.x, uy=(double)b.y-a.y;
    double vx=(double)c.x-a.x, vy=(double)c.y-a.y;
    double nu=hypot(ux,uy), nv=hypot(vx,vy);
    if (!(nu>1e-12 && nv>1e-12)) return 0.0;
    {
        double cs=(ux*vx+uy*vy)/(nu*nv);
        if (cs>1.0) cs=1.0;
        if (cs<-1.0) cs=-1.0;
        return acos(cs);
    }
}

static double energy(const geomlib_trimesh *m, const geomlib_vec2 *uv,
                     const double *target)
{
    size_t t;
    double e=0.0;
    for (t=0;t<m->triangle_count;++t) {
        uint32_t a=m->indices[t*3], b=m->indices[t*3+1], c=m->indices[t*3+2];
        double ea=target[t*3]-angle2(uv[a],uv[b],uv[c]);
        double eb=target[t*3+1]-angle2(uv[b],uv[c],uv[a]);
        double ec=target[t*3+2]-angle2(uv[c],uv[a],uv[b]);
        double ax=(double)uv[b].x-uv[a].x, ay=(double)uv[b].y-uv[a].y;
        double bx=(double)uv[c].x-uv[a].x, by=(double)uv[c].y-uv[a].y;
        double area=ax*by-ay*bx;
        e += ea*ea+eb*eb+ec*ec;
        if (area < 0.0) e += 100.0*area*area;
    }
    return e;
}

float geomlib_abf_angle_error(const geomlib_trimesh *mesh,
                              const geomlib_uv_coords *uvs)
{
    size_t t;
    double mx=0.0;
    if (!mesh || !uvs || !uvs->uvs) return INFINITY;
    for (t=0;t<mesh->triangle_count;++t) {
        uint32_t a=mesh->indices[t*3], b=mesh->indices[t*3+1], c=mesh->indices[t*3+2];
        geomlib_vec3 p=mesh->vertices[a], q=mesh->vertices[b], r=mesh->vertices[c];
        double d;
        if (a>=uvs->vertex_count||b>=uvs->vertex_count||c>=uvs->vertex_count)
            return INFINITY;
        d=fabs(angle3(p,q,r)-angle2(uvs->uvs[a],uvs->uvs[b],uvs->uvs[c]));
        if (d>mx) mx=d;
        d=fabs(angle3(q,r,p)-angle2(uvs->uvs[b],uvs->uvs[c],uvs->uvs[a]));
        if (d>mx) mx=d;
        d=fabs(angle3(r,p,q)-angle2(uvs->uvs[c],uvs->uvs[a],uvs->uvs[b]));
        if (d>mx) mx=d;
    }
    return (float)mx;
}

int geomlib_parameterize_abf(const geomlib_trimesh *mesh,
                             uint32_t pin0, uint32_t pin1,
                             const geomlib_abf_params *params,
                             geomlib_uv_coords *out)
{
    geomlib_abf_params p=params ? *params : geomlib_abf_params_default();
    double *target=NULL;
    size_t t, iter, v;
    int rc;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->triangle_count || pin0>=mesh->vertex_count ||
        pin1>=mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;
    *out=geomlib_uv_coords_create(mesh->vertex_count,mesh->triangle_count);
    if (!out->uvs) return GEOMLIB_MESH_ERROR_MEMORY;
    rc=geomlib_parameterize_lscm_free(mesh,pin0,pin1,out);
    if (rc!=GEOMLIB_MESH_SUCCESS) { geomlib_uv_coords_free(out); return rc; }
    if (p.iterations<1) p.iterations=1;
    if (p.step<=0.0f) p.step=0.05f;

    target=(double*)malloc(mesh->triangle_count*3*sizeof(double));
    if (!target) { geomlib_uv_coords_free(out); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (t=0;t<mesh->triangle_count;++t) {
        uint32_t a=mesh->indices[t*3], b=mesh->indices[t*3+1], c=mesh->indices[t*3+2];
        target[t*3]=angle3(mesh->vertices[a],mesh->vertices[b],mesh->vertices[c]);
        target[t*3+1]=angle3(mesh->vertices[b],mesh->vertices[c],mesh->vertices[a]);
        target[t*3+2]=angle3(mesh->vertices[c],mesh->vertices[a],mesh->vertices[b]);
    }

    for (iter=0;iter<(size_t)p.iterations;++iter) {
        if (geomlib_abf_angle_error(mesh,out) <= p.tolerance) break;
        for (v=0;v<mesh->vertex_count;++v) {
            geomlib_vec2 orig;
            double gu=0.0, gv=0.0, e0, eu;
            float h, alpha;
            int sub;
            if (v==pin0 || v==pin1 || v>=out->vertex_count) continue;
            orig=out->uvs[v];
            h=0.001f;
            e0=energy(mesh,out->uvs,target);
            out->uvs[v]=geomlib_vec2_create(orig.x+h,orig.y);
            eu=energy(mesh,out->uvs,target);
            out->uvs[v]=geomlib_vec2_create(orig.x-h,orig.y);
            gu=(eu-energy(mesh,out->uvs,target))/(2.0*h);
            out->uvs[v]=geomlib_vec2_create(orig.x,orig.y+h);
            eu=energy(mesh,out->uvs,target);
            out->uvs[v]=geomlib_vec2_create(orig.x,orig.y-h);
            gv=(eu-energy(mesh,out->uvs,target))/(2.0*h);
            out->uvs[v]=orig;
            if (!isfinite(gu) || !isfinite(gv)) continue;
            /* Backtracking: accept only an energy-decreasing update. */
            alpha=p.step;
            for (sub=0;sub<12;++sub) {
                out->uvs[v]=geomlib_vec2_create(
                    orig.x-(float)(alpha*gu), orig.y-(float)(alpha*gv));
                if (isfinite(out->uvs[v].x) && isfinite(out->uvs[v].y) &&
                    energy(mesh,out->uvs,target) < e0) break;
                alpha*=0.5f;
            }
            if (sub==12) out->uvs[v]=orig;
        }
    }
    free(target);
    return GEOMLIB_MESH_SUCCESS;
}
