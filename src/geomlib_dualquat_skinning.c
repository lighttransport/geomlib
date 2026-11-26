/*
 * geomlib - Dual-quaternion skinning (see geomlib_dualquat_skinning.h)
 *
 * Clean-room implementation of Kavan et al., "Skinning with Dual
 * Quaternions", SIGGRAPH 2007 (unit dual-quaternion blending).
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_dualquat_skinning.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

typedef struct { double x,y,z,w; } quat;

static quat q_mul(quat a, quat b)
{
    quat r;
    r.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
    r.x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
    r.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
    r.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
    return r;
}

static double q_norm(quat q)
{
    return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
}

/* Extract a unit dual quaternion (real part q0, dual part qe) from a
 * column-major rigid 4x4 transform. */
static void mat_to_dq(const geomlib_mat4 *m, quat *q0, quat *qe)
{
    const float *a = (const float *)m->m;
    double tr = a[0] + a[5] + a[10];
    double w, x, y, z, s;
    /* Standard quaternion extraction via the largest diagonal term. */
    {
        double m00=a[0], m11=a[5], m22=a[10];
        if (m00 > m11 && m00 > m22) {
            s = sqrt(1.0 + m00 - m11 - m22) * 2.0;
            w = (a[6] - a[9])/s;
            x = 0.25*s;
            y = (a[1] + a[4])/s;
            z = (a[8] + a[2])/s;
        } else if (m11 > m22) {
            s = sqrt(1.0 + m11 - m00 - m22) * 2.0;
            w = (a[8] - a[2])/s;
            x = (a[1] + a[4])/s;
            y = 0.25*s;
            z = (a[6] + a[9])/s;
        } else if (tr > 0.0) {
            s = sqrt(tr + 1.0)*2.0;
            w = 0.25*s;
            x = (a[6] - a[9])/s;
            y = (a[8] - a[2])/s;
            z = (a[1] - a[4])/s;
        } else {
            s = sqrt(1.0 + m22 - m00 - m11) * 2.0;
            w = (a[1] - a[4])/s;
            x = (a[8] + a[2])/s;
            y = (a[6] + a[9])/s;
            z = 0.25*s;
        }
    }
    q0->w=w; q0->x=x; q0->y=y; q0->z=z;
    {
        double n = q_norm(*q0);
        if (n > 1e-12) { q0->w/=n; q0->x/=n; q0->y/=n; q0->z/=n; }
        else { q0->w=1; q0->x=q0->y=q0->z=0; }
    }
    {
        quat tq = {a[12], a[13], a[14], 0.0};
        *qe = q_mul(tq, *q0);
        qe->x *= 0.5; qe->y *= 0.5; qe->z *= 0.5; qe->w *= 0.5;
    }
}

static geomlib_vec3 dq_transform(quat q, quat qe, geomlib_vec3 v)
{
    geomlib_vec3 cv, res, t;
    double vx=v.x, vy=v.y, vz=v.z;
    /* v + 2 * cross(q.xyz, cross(q.xyz, v) + q.w * v) */
    double cx = q.y*vz - q.z*vy;
    double cy = q.z*vx - q.x*vz;
    double cz = q.x*vy - q.y*vx;
    cx += q.w*vx; cy += q.w*vy; cz += q.w*vz;
    double rx = vx + 2.0*(q.y*cz - q.z*cy);
    double ry = vy + 2.0*(q.z*cx - q.x*cz);
    double rz = vz + 2.0*(q.x*cy - q.y*cx);
    res = geomlib_vec3_create((float)rx, (float)ry, (float)rz);
    /* translation = 2 * (qe * conjugate(q)).xyz */
    {
        quat cq = { -q.x, -q.y, -q.z, q.w };
        quat tq = q_mul(qe, cq);
        t = geomlib_vec3_create((float)(2.0*tq.x), (float)(2.0*tq.y), (float)(2.0*tq.z));
    }
    cv = geomlib_vec3_add(res, t);
    return cv;
}

int geomlib_skin_mesh_dualquat(const geomlib_trimesh *input,
                               const geomlib_skin_weight *weights,
                               const geomlib_mat4 *bone_transforms,
                               size_t bone_count,
                               geomlib_trimesh *output)
{
    geomlib_trimesh result = {0};
    size_t i, k;
    int rc;
    quat *q0, *qe;

    if (!input || !output || !input->vertices || !input->indices ||
        !weights || !bone_transforms || bone_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    q0 = (quat *)malloc(bone_count*sizeof(quat));
    qe = (quat *)malloc(bone_count*sizeof(quat));
    if (!q0 || !qe) { free(q0); free(qe); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (k = 0; k < bone_count; ++k) mat_to_dq(&bone_transforms[k], &q0[k], &qe[k]);

    rc = geomlib_trimesh_init(&result, input->vertex_count, input->triangle_count);
    if (rc != GEOMLIB_MESH_SUCCESS) { free(q0); free(qe); return rc; }

    for (i = 0; i < input->vertex_count; ++i) {
        quat b0 = {0,0,0,0}, be = {0,0,0,0};
        double wsum = 0.0;
        for (k = 0; k < GEOMLIB_SKIN_MAX_BONES; ++k) {
            float w = weights[i].weights[k];
            uint32_t b = weights[i].bone_indices[k];
            quat r0, re;
            if (!(w > 0.0f) || b >= bone_count) continue;
            r0 = q0[b]; re = qe[b];
            if (b0.w*r0.w + b0.x*r0.x + b0.y*r0.y + b0.z*r0.z < 0.0) {
                r0.w=-r0.w; r0.x=-r0.x; r0.y=-r0.y; r0.z=-r0.z;
                re.w=-re.w; re.x=-re.x; re.y=-re.y; re.z=-re.z;
            }
            b0.w += w*r0.w; b0.x += w*r0.x; b0.y += w*r0.y; b0.z += w*r0.z;
            be.w += w*re.w; be.x += w*re.x; be.y += w*re.y; be.z += w*re.z;
            wsum += w;
        }
        if (wsum > 1e-12) {
            double n = q_norm(b0);
            if (n > 1e-12) {
                b0.w/=n; b0.x/=n; b0.y/=n; b0.z/=n;
                be.w/=n; be.x/=n; be.y/=n; be.z/=n;
            }
            geomlib_trimesh_add_vertex(&result,
                dq_transform(b0, be, input->vertices[i]), NULL, NULL);
        } else {
            geomlib_trimesh_add_vertex(&result, input->vertices[i], NULL, NULL);
        }
    }
    for (i = 0; i < input->triangle_count; ++i)
        geomlib_trimesh_add_triangle(&result, input->indices[i*3],
                                     input->indices[i*3+1], input->indices[i*3+2]);
    *output = result;
    free(q0); free(qe);
    return GEOMLIB_MESH_SUCCESS;
}
