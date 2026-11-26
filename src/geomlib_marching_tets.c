/*
 * geomlib - Marching tetrahedra (see geomlib_marching_tets.h)
 *
 * Clean-room implementation. Each grid cube is decomposed into six
 * tetrahedra sharing the body diagonal corner 0->6; a tetrahedron with 1 or
 * 3 negative corners produces one triangle, and 2 negative corners produce
 * two triangles. Vertices at grid-edge crossings are welded afterwards.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_marching_tets.h"
#include "geomlib/geomlib_repair.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Cube corners use the same numbering as geomlib_sdf.c. */
static const int cube_corners[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
};

/* Six tetrahedra sharing the body diagonal 0-6. */
static const int cube_tets[6][4] = {
    {0,5,1,6}, {0,1,2,6}, {0,2,3,6}, {0,3,7,6}, {0,7,4,6}, {0,4,5,6}
};

typedef struct {
    geomlib_vec3 *v;
    size_t count, cap;
    uint32_t *idx;
    size_t icount, icap;
} soup;

static void add_vertex(soup *s, geomlib_vec3 p, uint32_t *out)
{
    if (s->count == s->cap) {
        size_t nc = s->cap ? s->cap*2 : 4096;
        geomlib_vec3 *nv = (geomlib_vec3 *)realloc(s->v, nc*sizeof(geomlib_vec3));
        if (!nv) return;
        s->v = nv;
        s->cap = nc;
    }
    s->v[s->count] = p;
    *out = (uint32_t)s->count++;
}

/* Emit a triangle whose normal points from the inside corners toward the
 * outside corners of its tetrahedron (outward for a negative-inside field). */
static void add_triangle_oriented(soup *s, uint32_t a, uint32_t b, uint32_t c,
                                  geomlib_vec3 inside_center, geomlib_vec3 outside_center)
{
    geomlib_vec3 n = geomlib_vec3_cross(
        geomlib_vec3_sub(s->v[b], s->v[a]), geomlib_vec3_sub(s->v[c], s->v[a]));
    geomlib_vec3 d = geomlib_vec3_sub(outside_center, inside_center);
    if (geomlib_vec3_dot(n, d) < 0.0f) {
        uint32_t t = b; b = c; c = t;
    }
    if (s->icount + 3 > s->icap) {
        size_t nc = s->icap ? s->icap*2 : 4096;
        uint32_t *ni = (uint32_t *)realloc(s->idx, nc*sizeof(uint32_t));
        if (!ni) return;
        s->idx = ni;
        s->icap = nc;
    }
    s->idx[s->icount++] = a;
    s->idx[s->icount++] = b;
    s->idx[s->icount++] = c;
}

static geomlib_vec3 corner_pos(const geomlib_sdf *sdf, size_t x, size_t y, size_t z,
                               int corner)
{
    geomlib_vec3 p;
    p.x = sdf->origin.x + ((float)x + cube_corners[corner][0]) * sdf->voxel_size;
    p.y = sdf->origin.y + ((float)y + cube_corners[corner][1]) * sdf->voxel_size;
    p.z = sdf->origin.z + ((float)z + cube_corners[corner][2]) * sdf->voxel_size;
    return p;
}

static geomlib_vec3 lerp_pos(geomlib_vec3 a, geomlib_vec3 b, double va, double vb)
{
    double t = (0.0 - va) / (vb - va);
    if (!(t > 0.0 && t < 1.0)) t = 0.5;
    return geomlib_vec3_add(geomlib_vec3_mul(a, (float)(1.0-t)),
                            geomlib_vec3_mul(b, (float)t));
}

int geomlib_sdf_marching_tets(const geomlib_sdf *sdf, geomlib_trimesh *out)
{
    soup s;
    size_t x, y, z, t;
    uint32_t res[3];

    if (!sdf || !out || !sdf->values) return GEOMLIB_MESH_ERROR_INVALID;
    res[0] = sdf->resolution[0];
    res[1] = sdf->resolution[1];
    res[2] = sdf->resolution[2];
    if (res[0] < 2 || res[1] < 2 || res[2] < 2) return GEOMLIB_MESH_ERROR_INVALID;

    memset(&s, 0, sizeof(s));
    for (z = 0; z + 1 < res[2]; ++z) {
        for (y = 0; y + 1 < res[1]; ++y) {
            for (x = 0; x + 1 < res[0]; ++x) {
                double v[8];
                int c;
                for (c = 0; c < 8; ++c) {
                    size_t cx = x + cube_corners[c][0];
                    size_t cy = y + cube_corners[c][1];
                    size_t cz = z + cube_corners[c][2];
                    v[c] = sdf->values[cx + res[0]*(cy + res[1]*cz)];
                }
                for (t = 0; t < 6; ++t) {
                    const int *tet = cube_tets[t];
                    int inside[4], nin = 0;
                    int a, b;
                    for (a = 0; a < 4; ++a) {
                        inside[a] = v[tet[a]] < 0.0;
                        if (inside[a]) nin++;
                    }
                    if (nin == 0 || nin == 4) continue;
                    {
                        int pairs[6][2], np = 0;
                        uint32_t ids[6];
                        for (a = 0; a < 4; ++a) {
                            for (b = a+1; b < 4; ++b) {
                                if (inside[a] == inside[b]) continue;
                                pairs[np][0] = a;
                                pairs[np][1] = b;
                                np++;
                            }
                        }
                        if (np != 3 && np != 4) continue;
                        for (a = 0; a < np; ++a) {
                            geomlib_vec3 pa = corner_pos(sdf, x,y,z,tet[pairs[a][0]]);
                            geomlib_vec3 pb = corner_pos(sdf, x,y,z,tet[pairs[a][1]]);
                            add_vertex(&s, lerp_pos(pa, pb,
                                       v[tet[pairs[a][0]]], v[tet[pairs[a][1]]]),
                                       &ids[a]);
                        }
                        {
                            geomlib_vec3 ic = {0,0,0}, oc = {0,0,0};
                            int ni = 0, no = 0;
                            for (a = 0; a < 4; ++a) {
                                geomlib_vec3 cp = corner_pos(sdf, x,y,z,tet[a]);
                                if (inside[a]) { ic = geomlib_vec3_add(ic, cp); ni++; }
                                else { oc = geomlib_vec3_add(oc, cp); no++; }
                            }
                            if (ni) ic = geomlib_vec3_mul(ic, 1.0f/(float)ni);
                            if (no) oc = geomlib_vec3_mul(oc, 1.0f/(float)no);
                            if (np == 3) {
                                add_triangle_oriented(&s, ids[0], ids[1], ids[2], ic, oc);
                            } else {
                                /* Order the four crossings into a cycle: two
                                 * crossing edges are neighbours when they
                                 * share a tetrahedron corner. */
                                uint32_t ordered[4];
                                int used[4] = {0,0,0,0};
                                int prev0 = pairs[0][0], prev1 = pairs[0][1];
                                ordered[0] = ids[0];
                                used[0] = 1;
                                for (a = 1; a < 4; ++a) {
                                    int found = -1;
                                    for (b = 0; b < 4; ++b) {
                                        if (used[b]) continue;
                                        if (pairs[b][0]==prev0 || pairs[b][0]==prev1 ||
                                            pairs[b][1]==prev0 || pairs[b][1]==prev1) {
                                            found = b;
                                            break;
                                        }
                                    }
                                    if (found < 0) found = a;
                                    ordered[a] = ids[found];
                                    used[found] = 1;
                                    prev0 = pairs[found][0];
                                    prev1 = pairs[found][1];
                                }
                                add_triangle_oriented(&s, ordered[0], ordered[1],
                                                      ordered[2], ic, oc);
                                add_triangle_oriented(&s, ordered[0], ordered[2],
                                                      ordered[3], ic, oc);
                            }
                        }
                    }
                }
            }
        }
    }

    if (s.icount == 0) {
        geomlib_trimesh_init(out, 0, 0);
        free(s.v); free(s.idx);
        return GEOMLIB_MESH_SUCCESS;
    }
    {
        int rc = geomlib_soup_to_mesh(s.v, s.count, s.idx, s.icount,
                                      sdf->voxel_size * 1e-3f, out);
        free(s.v); free(s.idx);
        return rc;
    }
}
