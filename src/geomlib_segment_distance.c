/*
 * geomlib - Segment-segment distance (see geomlib_segment_distance.h)
 *
 * Clean-room clamped closest-point iteration with exact parallel handling.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_segment_distance.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static float clamp01(float v)
{
    if (v<0.0f) return 0.0f;
    if (v>1.0f) return 1.0f;
    return v;
}

int geomlib_segment_closest_points(geomlib_vec3 a0, geomlib_vec3 a1,
                                   geomlib_vec3 b0, geomlib_vec3 b1,
                                   geomlib_segment_distance *out)
{
    geomlib_vec3 u=geomlib_vec3_sub(a1,a0),v=geomlib_vec3_sub(b1,b0);
    geomlib_vec3 w=geomlib_vec3_sub(a0,b0);
    float a=geomlib_vec3_dot(u,u),b=geomlib_vec3_dot(u,v),c=geomlib_vec3_dot(v,v);
    float d=geomlib_vec3_dot(u,w),e=geomlib_vec3_dot(v,w);
    float den=a*c-b*b,s,t;
    geomlib_vec3 pa,pb;

    if (!out) return GEOMLIB_MESH_ERROR_INVALID;
    if (den>1e-12f) {
        s=clamp01((b*e-c*d)/den);
        t=(b*s+e)/c;
        if (t<0.0f) { t=0.0f; s=clamp01(-d/a); }
        else if (t>1.0f) { t=1.0f; s=clamp01((b-d)/a); }
    } else {
        /* Parallel (or degenerate): try every endpoint-to-segment projection
         * and keep the closest pair. This handles collinear segments whose
         * ranges overlap or are disjoint. */
        float best=FLT_MAX,bs=0.0f,bt=0.0f;
        int k;
        for (k=0;k<2;++k) {
            float sk=(float)k;
            geomlib_vec3 pk=geomlib_vec3_add(a0,geomlib_vec3_mul(u,sk));
            float tt = c>1e-12f ? clamp01(geomlib_vec3_dot(
                geomlib_vec3_sub(pk,b0),v)/c) : 0.0f;
            geomlib_vec3 qb=geomlib_vec3_add(b0,geomlib_vec3_mul(v,tt));
            float dd=geomlib_vec3_length_squared(geomlib_vec3_sub(pk,qb));
            if (dd<best) { best=dd; bs=sk; bt=tt; }

            float tk=(float)k;
            geomlib_vec3 qk=geomlib_vec3_add(b0,geomlib_vec3_mul(v,tk));
            float ss = a>1e-12f ? clamp01(geomlib_vec3_dot(
                geomlib_vec3_sub(qk,a0),u)/a) : 0.0f;
            geomlib_vec3 pa2=geomlib_vec3_add(a0,geomlib_vec3_mul(u,ss));
            float dd2=geomlib_vec3_length_squared(geomlib_vec3_sub(qk,pa2));
            if (dd2<best) { best=dd2; bs=ss; bt=tk; }
        }
        s=bs; t=bt;
    }
    if (a<1e-12f) s=0.0f;
    pa=geomlib_vec3_add(a0,geomlib_vec3_mul(u,s));
    pb=geomlib_vec3_add(b0,geomlib_vec3_mul(v,t));
    out->point_a=pa;
    out->point_b=pb;
    out->s=s; out->t=t;
    out->distance=geomlib_vec3_length(geomlib_vec3_sub(pa,pb));
    return GEOMLIB_MESH_SUCCESS;
}
