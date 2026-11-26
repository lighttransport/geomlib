/*
 * geomlib - Polygon triangulation with holes (see geomlib_polytri_holes.h)
 *
 * Clean-room implementation: build a constrained Delaunay triangulation of
 * all boundary vertices with every outer and hole edge constrained, then drop
 * triangles whose centroid lies inside a hole.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polytri_holes.h"
#include "geomlib/geomlib_cdt2d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_polytri_result geomlib_polytri_result_create(void)
{
    geomlib_polytri_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_polytri_result_free(geomlib_polytri_result *r)
{
    if (!r) return;
    free(r->triangles);
    free(r->vertices);
    memset(r, 0, sizeof(*r));
}

static int inside_polygon(geomlib_vec2 q, const geomlib_vec2 *p, size_t n)
{
    size_t i;
    int inside = 0;
    for (i = 0; i < n; ++i) {
        size_t j=(i+1)%n;
        if (((p[i].y>q.y)!=(p[j].y>q.y)) &&
            ((double)q.x < ((double)p[j].x-p[i].x)*((double)q.y-p[i].y)/
                           ((double)p[j].y-p[i].y)+p[i].x))
            inside=!inside;
    }
    return inside;
}

int geomlib_triangulate_with_holes(const geomlib_vec2 *outer, size_t no,
                                   const geomlib_vec2 *const *holes,
                                   const size_t *hole_sizes, size_t nh,
                                   geomlib_polytri_result *out)
{
    geomlib_polytri_result result = geomlib_polytri_result_create();
    geomlib_cdt2d cdt = geomlib_cdt2d_create();
    geomlib_vec2 *verts = NULL;
    uint32_t *cons = NULL;
    size_t nv = 0, nc = 0, i, h;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !outer || no < 3) return rc;
    if (nh && (!holes || !hole_sizes)) return rc;

    for (h = 0; h <= nh; ++h)
        if (h == 0) nv += no; else nv += hole_sizes[h-1];
    verts = (geomlib_vec2 *)malloc(nv*sizeof(geomlib_vec2));
    cons = (uint32_t *)malloc(2*nv*sizeof(uint32_t));
    if (!verts || !cons) goto cleanup;

    nv = 0;
    for (i = 0; i < no; ++i) verts[nv++] = outer[i];
    for (h = 0; h < nh; ++h)
        for (i = 0; i < hole_sizes[h]; ++i) verts[nv++] = holes[h][i];

    nv = 0;
    for (i = 0; i < no; ++i) {
        cons[nc++] = (uint32_t)i;
        cons[nc++] = (uint32_t)((i+1)%no);
    }
    nv = no;
    for (h = 0; h < nh; ++h) {
        for (i = 0; i < hole_sizes[h]; ++i) {
            cons[nc++] = (uint32_t)(nv+i);
            cons[nc++] = (uint32_t)(nv+(i+1)%hole_sizes[h]);
        }
        nv += hole_sizes[h];
    }
    nv = 0;
    for (h = 0; h <= nh; ++h)
        if (h == 0) nv += no; else nv += hole_sizes[h-1];

    rc = geomlib_cdt2d_build(verts, nv, cons, nc/2, &cdt);
    if (rc != GEOMLIB_MESH_SUCCESS) goto cleanup;

    result.triangles = (uint32_t *)malloc(cdt.triangle_count*3*sizeof(uint32_t));
    if (!result.triangles) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
    for (i = 0; i < cdt.triangle_count; ++i) {
        uint32_t a=cdt.triangles[i*3], b=cdt.triangles[i*3+1], c=cdt.triangles[i*3+2];
        geomlib_vec2 mid = geomlib_vec2_create(
            (verts[a].x+verts[b].x+verts[c].x)/3.0f,
            (verts[a].y+verts[b].y+verts[c].y)/3.0f);
        int in_hole = 0;
        size_t off = no, q;
        for (q = 0; q < nh; ++q) {
            if (inside_polygon(mid, &verts[off], hole_sizes[q])) { in_hole=1; break; }
            off += hole_sizes[q];
        }
        if (in_hole) continue;
        result.triangles[result.triangle_count*3]=a;
        result.triangles[result.triangle_count*3+1]=b;
        result.triangles[result.triangle_count*3+2]=c;
        result.triangle_count++;
    }
    result.vertices = verts; verts=NULL;
    result.vertex_count = nv;
    geomlib_cdt2d_free(&cdt);
    free(cons);
    *out = result;
    return GEOMLIB_MESH_SUCCESS;

cleanup:
    free(verts); free(cons);
    geomlib_cdt2d_free(&cdt);
    geomlib_polytri_result_free(&result);
    return rc;
}
