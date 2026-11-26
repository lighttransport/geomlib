/*
 * geomlib - Mesh contour lines (see geomlib_mesh_contour.h)
 *
 * Clean-room marching-triangles on a triangle mesh. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_mesh_contour.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_mesh_contour geomlib_mesh_contour_create(void)
{
    geomlib_mesh_contour c;
    memset(&c, 0, sizeof(c));
    return c;
}

void geomlib_mesh_contour_free(geomlib_mesh_contour *c)
{
    if (!c) return;
    free(c->segments);
    memset(c, 0, sizeof(*c));
}

static int add_segment(geomlib_mesh_contour *c, geomlib_vec3 a, geomlib_vec3 b,
                       size_t *cap)
{
    if (geomlib_vec3_length(geomlib_vec3_sub(a,b)) < 1e-9f) return 0;
    if (c->segment_count==*cap) {
        size_t nc=*cap?(*cap)*2:64;
        geomlib_vec3 *ns=(geomlib_vec3*)realloc(c->segments,nc*2*sizeof(geomlib_vec3));
        if (!ns) return -1;
        c->segments=ns; *cap=nc;
    }
    c->segments[c->segment_count*2]=a;
    c->segments[c->segment_count*2+1]=b;
    c->segment_count++;
    return 0;
}

int geomlib_mesh_contour_lines(const geomlib_trimesh *mesh,
                               const float *values, float iso,
                               geomlib_mesh_contour *out)
{
    geomlib_mesh_contour r=geomlib_mesh_contour_create();
    size_t cap=0, t, i;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !values || !mesh->triangle_count) return rc;
    for (i=0;i<mesh->vertex_count;++i)
        if (!isfinite(values[i])) return rc;

    for (t=0;t<mesh->triangle_count;++t) {
        uint32_t v[3]={mesh->indices[t*3],mesh->indices[t*3+1],mesh->indices[t*3+2]};
        geomlib_vec3 p[3];
        float f[3];
        geomlib_vec3 hits[4];
        size_t nh=0;
        int k;
        for (k=0;k<3;++k) {
            if (v[k]>=mesh->vertex_count) goto cleanup;
            p[k]=mesh->vertices[v[k]];
            f[k]=values[v[k]];
        }
        for (k=0;k<3;++k) {
            uint32_t a=(uint32_t)k,b=(uint32_t)((k+1)%3);
            float fa=f[a]-iso, fb=f[b]-iso;
            if (fabs(fa)<1e-9f) {
                if (nh<4) hits[nh++]=p[a];
                continue;
            }
            if (fa*fb<0.0f) {
                float u=fa/(fa-fb);
                if (nh<4) hits[nh++]=geomlib_vec3_add(
                    p[a],geomlib_vec3_mul(geomlib_vec3_sub(p[b],p[a]),u));
            }
        }
        if (nh==2) {
            if (add_segment(&r,hits[0],hits[1],&cap)<0) goto cleanup;
        } else if (nh==4) {
            /* Two ambiguous crossings: connect by proximity to avoid
             * crossing diagonals. */
            double d01=geomlib_vec3_length_squared(geomlib_vec3_sub(hits[0],hits[1]));
            double d02=geomlib_vec3_length_squared(geomlib_vec3_sub(hits[0],hits[2]));
            double d03=geomlib_vec3_length_squared(geomlib_vec3_sub(hits[0],hits[3]));
            if (d01<=d02 && d01<=d03) {
                if (add_segment(&r,hits[0],hits[1],&cap)<0) goto cleanup;
                if (add_segment(&r,hits[2],hits[3],&cap)<0) goto cleanup;
            } else if (d02<=d03) {
                if (add_segment(&r,hits[0],hits[2],&cap)<0) goto cleanup;
                if (add_segment(&r,hits[1],hits[3],&cap)<0) goto cleanup;
            } else {
                if (add_segment(&r,hits[0],hits[3],&cap)<0) goto cleanup;
                if (add_segment(&r,hits[1],hits[2],&cap)<0) goto cleanup;
            }
        }
    }
    *out=r;
    return GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_mesh_contour_free(&r);
    return rc;
}
