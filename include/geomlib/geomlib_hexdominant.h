/*
 * geomlib - Hex-dominant meshing via tetrahedral cell merging
 *
 * Detects groups of five tetrahedra sharing a body diagonal (the standard
 * 5-tet split of a hexahedron) and merges each group into one hex cell.
 */

#ifndef GEOMLIB_HEXDOMINANT_H
#define GEOMLIB_HEXDOMINANT_H

#include "geomlib_tetramesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *hexes;        /* 8 indices per hex cell */
    size_t hex_count;
    uint32_t *tets;         /* 4 indices per remaining tetrahedron */
    size_t tet_count;
    geomlib_vec3 *vertices; /* shared with the input tetrahedral mesh */
    size_t vertex_count;
} geomlib_hexmesh;

geomlib_hexmesh geomlib_hexmesh_create(void);
void geomlib_hexmesh_free(geomlib_hexmesh *m);

/* Convert a tetrahedral mesh to a mixed hex/tet mesh. Values are borrowed
 * from `input`; free only the index arrays with geomlib_hexmesh_free. */
int geomlib_hex_dominant(const geomlib_tetramesh *input,
                         geomlib_hexmesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_HEXDOMINANT_H */
