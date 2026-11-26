/*
 * geomlib - Tetrahedral cage deformation (see geomlib_cage3d.h)
 *
 * Clean-room volume-ratio barycentric coordinates. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_cage3d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double tet_volume(geomlib_vec3 a, geomlib_vec3 b,
                         geomlib_vec3 c, geomlib_vec3 d)
{
    geomlib_vec3 u=geomlib_vec3_sub(b,a),v=geomlib_vec3_sub(c,a),
                 w=geomlib_vec3_sub(d,a);
    return geomlib_vec3_dot(u,geomlib_vec3_cross(v,w))/6.0;
}

int geomlib_tetra_weights(const geomlib_vec3 cage[4], geomlib_vec3 q,
                          double weights[4])
{
    double total;
    if (!cage || !weights) return GEOMLIB_MESH_ERROR_INVALID;
    total=tet_volume(cage[0],cage[1],cage[2],cage[3]);
    if (fabs(total)<1e-14) return GEOMLIB_MESH_ERROR_INVALID;
    weights[0]=tet_volume(q,cage[1],cage[2],cage[3])/total;
    weights[1]=tet_volume(cage[0],q,cage[2],cage[3])/total;
    weights[2]=tet_volume(cage[0],cage[1],q,cage[3])/total;
    weights[3]=tet_volume(cage[0],cage[1],cage[2],q)/total;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_deform_tetra_cage(const geomlib_vec3 cage[4],
                              const geomlib_vec3 deformed[4],
                              const geomlib_vec3 *points, size_t count,
                              geomlib_vec3 *out)
{
    double w[4];
    size_t i;
    int k,rc;
    if (!cage||!deformed||!points||!out) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<count;++i) {
        rc=geomlib_tetra_weights(cage,points[i],w);
        if (rc!=GEOMLIB_MESH_SUCCESS) return rc;
        out[i]=geomlib_vec3_create(0,0,0);
        for (k=0;k<4;++k)
            out[i]=geomlib_vec3_add(out[i],geomlib_vec3_mul(deformed[k],(float)w[k]));
    }
    return GEOMLIB_MESH_SUCCESS;
}
