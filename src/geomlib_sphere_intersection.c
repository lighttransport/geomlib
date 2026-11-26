/*
 * geomlib - Sphere intersection (see geomlib_sphere_intersection.h)
 *
 * Clean-room analytic solution. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sphere_intersection.h"

#include <float.h>
#include <math.h>
#include <string.h>

geomlib_sphere_intersection geomlib_sphere_sphere_intersection(
    geomlib_vec3 c0, float r0, geomlib_vec3 c1, float r1)
{
    geomlib_sphere_intersection r;
    geomlib_vec3 d=geomlib_vec3_sub(c1,c0);
    float dist=geomlib_vec3_length(d),a;

    memset(&r,0,sizeof(r));
    r.kind=GEOMLIB_SPHERE_NONE;
    r.radius=0.0f;
    if (!(r0>=0.0f && r1>=0.0f)) return r;
    if (dist<1e-9f) {
        if (fabsf(r0-r1)<1e-6f) { r.kind=GEOMLIB_SPHERE_CONTAINED; r.center=c0; }
        return r;
    }
    r.normal=geomlib_vec3_mul(d,1.0f/dist);
    if (dist>r0+r1+1e-6f || dist<fabsf(r0-r1)-1e-6f) return r;
    a=(r0*r0-r1*r1+dist*dist)/(2.0f*dist);
    r.center=geomlib_vec3_add(c0,geomlib_vec3_mul(r.normal,a));
    {
        double rr=r0*r0-(double)a*a;
        if (rr<-1e-6f) return r;
        r.radius=(float)sqrt(fmax(0.0,rr));
    }
    r.kind=(r.radius<1e-5f) ? GEOMLIB_SPHERE_POINT : GEOMLIB_SPHERE_CIRCLE;
    return r;
}
