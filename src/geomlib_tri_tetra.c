/*
 * geomlib - Triangle / tetrahedron intersection (see geomlib_tri_tetra.h)
 *
 * Clean-room closest-feature enumeration. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_tri_tetra.h"

#include <stddef.h>
#include <math.h>
#include "geomlib/geomlib_tetramesh.h"
#include "geomlib/geomlib_segment_triangle.h"

static void tetra_face(const geomlib_vec3 v[4], int f,
                       geomlib_vec3 *a, geomlib_vec3 *b, geomlib_vec3 *c)
{
    static const int idx[4][3] = {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
    *a = v[idx[f][0]];
    *b = v[idx[f][1]];
    *c = v[idx[f][2]];
}

static int point_in_tri3d(geomlib_vec3 p, geomlib_vec3 a,
                          geomlib_vec3 b, geomlib_vec3 c)
{
    geomlib_vec3 e1 = geomlib_vec3_sub(b, a), e2 = geomlib_vec3_sub(c, a);
    geomlib_vec3 n = geomlib_vec3_cross(e1, e2);
    geomlib_vec3 q = geomlib_vec3_sub(p, a);
    float d00 = geomlib_vec3_dot(e1, e1);
    float d01 = geomlib_vec3_dot(e1, e2);
    float d11 = geomlib_vec3_dot(e2, e2);
    float den = d00 * d11 - d01 * d01;
    float u, v;
    if (!(den > 1e-14f)) return 0;
    if (fabsf(geomlib_vec3_dot(n, q)) > 1e-5f *
        (geomlib_vec3_length(n) + geomlib_vec3_length(q) + 1.0f)) return 0;
    u = (d11 * geomlib_vec3_dot(e1, q) - d01 * geomlib_vec3_dot(e2, q)) / den;
    v = (d00 * geomlib_vec3_dot(e2, q) - d01 * geomlib_vec3_dot(e1, q)) / den;
    return u >= -1e-6f && v >= -1e-6f && u + v <= 1.0f + 1e-6f;
}

int geomlib_triangle_tetrahedron_intersect(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2, geomlib_vec3 v3)
{
    geomlib_vec3 tri[3], tet[4], f0, f1, f2;
    int i, f;

    tri[0] = a0; tri[1] = a1; tri[2] = a2;
    tet[0] = v0; tet[1] = v1; tet[2] = v2; tet[3] = v3;

    for (i = 0; i < 3; ++i)
        if (geomlib_point_in_tetrahedron(tri[i], v0, v1, v2, v3)) return 1;
    for (i = 0; i < 4; ++i)
        if (point_in_tri3d(tet[i], a0, a1, a2)) return 1;

    for (f = 0; f < 4; ++f) {
        tetra_face(tet, f, &f0, &f1, &f2);
        for (i = 0; i < 3; ++i)
            if (geomlib_segment_triangle_intersect(
                    tri[i], tri[(i + 1) % 3], f0, f1, f2, NULL)) return 1;
    }
    {
        static const int edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        for (i = 0; i < 6; ++i)
            if (geomlib_segment_triangle_intersect(
                    tet[edges[i][0]], tet[edges[i][1]], a0, a1, a2, NULL))
                return 1;
    }
    return 0;
}
