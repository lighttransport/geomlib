/*
 * geomlib - Tetrahedron / tetrahedron intersection (see geomlib_tetra_tetra.h)
 *
 * Clean-room closest-feature enumeration. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_tetra_tetra.h"
#include "geomlib/geomlib_tetramesh.h"
#include "geomlib/geomlib_segment_triangle.h"

#include <stddef.h>

static const int TETRA_FACES[4][3] = {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
static const int TETRA_EDGES[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

static int edges_against_faces(const geomlib_vec3 t[4], const geomlib_vec3 o[4])
{
    int i, f;
    for (i = 0; i < 6; ++i)
        for (f = 0; f < 4; ++f)
            if (geomlib_segment_triangle_intersect(
                    t[TETRA_EDGES[i][0]], t[TETRA_EDGES[i][1]],
                    o[TETRA_FACES[f][0]], o[TETRA_FACES[f][1]],
                    o[TETRA_FACES[f][2]], NULL)) return 1;
    return 0;
}

int geomlib_tetrahedra_intersect(
    geomlib_vec3 a0, geomlib_vec3 a1, geomlib_vec3 a2, geomlib_vec3 a3,
    geomlib_vec3 b0, geomlib_vec3 b1, geomlib_vec3 b2, geomlib_vec3 b3)
{
    geomlib_vec3 a[4], b[4];
    int i;

    a[0] = a0; a[1] = a1; a[2] = a2; a[3] = a3;
    b[0] = b0; b[1] = b1; b[2] = b2; b[3] = b3;

    for (i = 0; i < 4; ++i) {
        if (geomlib_point_in_tetrahedron(a[i], b0, b1, b2, b3)) return 1;
        if (geomlib_point_in_tetrahedron(b[i], a0, a1, a2, a3)) return 1;
    }
    if (edges_against_faces(a, b)) return 1;
    if (edges_against_faces(b, a)) return 1;
    return 0;
}
