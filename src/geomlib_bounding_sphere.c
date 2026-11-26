/*
 * geomlib - Ritter bounding sphere (see geomlib_bounding_sphere.h)
 *
 * Clean-room implementation: pick the farthest pair of extreme points along
 * the three axes, seed with that diameter, then grow the sphere point by point.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_bounding_sphere.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_bounding_sphere_points(const geomlib_vec3 *points, size_t count,
                                   geomlib_bounding_sphere *out)
{
    size_t i, p=0;
    geomlib_vec3 mn[3], mx[3], d;
    float best=-1.0f, rad;

    if (!out || !points || !count) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<3;++i) { mn[i]=points[0]; mx[i]=points[0]; }
    for (i=1;i<count;++i) {
        if (points[i].x<mn[0].x) mn[0]=points[i];
        if (points[i].x>mx[0].x) mx[0]=points[i];
        if (points[i].y<mn[1].y) mn[1]=points[i];
        if (points[i].y>mx[1].y) mx[1]=points[i];
        if (points[i].z<mn[2].z) mn[2]=points[i];
        if (points[i].z>mx[2].z) mx[2]=points[i];
    }
    for (i=0;i<3;++i) {
        d=geomlib_vec3_sub(mx[i],mn[i]);
        float len=geomlib_vec3_length(d);
        if (len>best) { best=len; p=i; }
    }
    out->center=geomlib_vec3_mul(geomlib_vec3_add(mn[p],mx[p]),0.5f);
    for (i=0;i<count;++i) {
        float dist=geomlib_vec3_length(geomlib_vec3_sub(points[i],out->center));
        if (dist>best) { best=dist; p=i; }
    }
    out->center=points[p];
    rad=0.0f;
    for (i=0;i<count;++i) {
        float dist=geomlib_vec3_length(geomlib_vec3_sub(points[i],out->center));
        if (dist>rad) { rad=dist; p=i; }
    }
    if (rad<=1e-12f) { out->radius=0.0f; return GEOMLIB_MESH_SUCCESS; }
    out->center=geomlib_vec3_mul(geomlib_vec3_add(out->center,points[p]),0.5f);
    rad=geomlib_vec3_length(geomlib_vec3_sub(points[p],out->center));
    out->radius=rad;

    /* Grow the sphere to enclose every point. */
    for (i=0;i<count;++i) {
        geomlib_vec3 off=geomlib_vec3_sub(points[i],out->center);
        float dist=geomlib_vec3_length(off);
        if (dist<=out->radius+1e-6f) continue;
        {
            float nr=0.5f*(out->radius+dist);
            float t=(dist-out->radius)/(2.0f*dist);
            out->center=geomlib_vec3_add(out->center,geomlib_vec3_mul(off,t));
            out->radius=nr;
        }
    }
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_bounding_sphere_mesh(const geomlib_trimesh *mesh,
                                 geomlib_bounding_sphere *out)
{
    if (!out || !mesh || !mesh->vertices || !mesh->vertex_count)
        return GEOMLIB_MESH_ERROR_INVALID;
    return geomlib_bounding_sphere_points(mesh->vertices,mesh->vertex_count,out);
}
