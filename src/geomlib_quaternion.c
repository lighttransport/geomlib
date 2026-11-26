/*
 * geomlib - Quaternions (see geomlib_quaternion.h)
 *
 * Clean-room unit-quaternion algebra and slerp. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_quaternion.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_quat geomlib_quat_identity(void)
{
    geomlib_quat q={1.0,0.0,0.0,0.0};
    return q;
}

geomlib_quat geomlib_quat_from_axis_angle(geomlib_vec3 axis, double angle)
{
    geomlib_quat q;
    float l=geomlib_vec3_length(axis);
    double s,c,h;
    if (!(l>1e-12f)) return geomlib_quat_identity();
    h=0.5*angle; s=sin(h); c=cos(h);
    q.w=c;
    q.x=s*axis.x/l; q.y=s*axis.y/l; q.z=s*axis.z/l;
    return q;
}

geomlib_quat geomlib_quat_mul(geomlib_quat a, geomlib_quat b)
{
    geomlib_quat q;
    q.w=a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z;
    q.x=a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y;
    q.y=a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x;
    q.z=a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w;
    return q;
}

geomlib_quat geomlib_quat_conjugate(geomlib_quat q)
{
    q.x=-q.x; q.y=-q.y; q.z=-q.z;
    return q;
}

geomlib_quat geomlib_quat_normalize(geomlib_quat q)
{
    double n=sqrt(q.w*q.w+q.x*q.x+q.y*q.y+q.z*q.z);
    if (!(n>1e-15)) return geomlib_quat_identity();
    q.w/=n; q.x/=n; q.y/=n; q.z/=n;
    return q;
}

double geomlib_quat_dot(geomlib_quat a, geomlib_quat b)
{
    return a.w*b.w+a.x*b.x+a.y*b.y+a.z*b.z;
}

geomlib_mat4 geomlib_quat_to_mat4(geomlib_quat q)
{
    geomlib_mat4 m=geomlib_mat4_identity();
    double w,x,y,z;
    q=geomlib_quat_normalize(q);
    w=q.w; x=q.x; y=q.y; z=q.z;
    m.m[0]=(float)(1.0-2.0*(y*y+z*z));
    m.m[1]=(float)(2.0*(x*y+w*z));
    m.m[2]=(float)(2.0*(x*z-w*y));
    m.m[4]=(float)(2.0*(x*y-w*z));
    m.m[5]=(float)(1.0-2.0*(x*x+z*z));
    m.m[6]=(float)(2.0*(y*z+w*x));
    m.m[8]=(float)(2.0*(x*z+w*y));
    m.m[9]=(float)(2.0*(y*z-w*x));
    m.m[10]=(float)(1.0-2.0*(x*x+y*y));
    return m;
}

geomlib_quat geomlib_quat_slerp(geomlib_quat a, geomlib_quat b, double t)
{
    double cosom,som,s0,s1,w;
    a=geomlib_quat_normalize(a);
    b=geomlib_quat_normalize(b);
    cosom=geomlib_quat_dot(a,b);
    if (cosom<0.0) { b.w=-b.w; b.x=-b.x; b.y=-b.y; b.z=-b.z; cosom=-cosom; }
    if (cosom>0.9995) {
        geomlib_quat q;
        q.w=a.w+t*(b.w-a.w);
        q.x=a.x+t*(b.x-a.x);
        q.y=a.y+t*(b.y-a.y);
        q.z=a.z+t*(b.z-a.z);
        return geomlib_quat_normalize(q);
    }
    w=acos(fmin(1.0,fmax(-1.0,cosom)));
    som=sin(w);
    s0=sin((1.0-t)*w)/som;
    s1=sin(t*w)/som;
    {
        geomlib_quat q;
        q.w=s0*a.w+s1*b.w;
        q.x=s0*a.x+s1*b.x;
        q.y=s0*a.y+s1*b.y;
        q.z=s0*a.z+s1*b.z;
        return q;
    }
}

geomlib_vec3 geomlib_quat_rotate(geomlib_quat q, geomlib_vec3 v)
{
    geomlib_mat4 m=geomlib_quat_to_mat4(q);
    geomlib_vec4 r=geomlib_mat4_mul_vec4(m,
        geomlib_vec4_create(v.x,v.y,v.z,1.0f));
    return geomlib_vec3_create(r.x,r.y,r.z);
}
