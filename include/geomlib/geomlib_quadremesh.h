/*
 * geomlib - Quad remeshing via triangle pairing
 *
 * Greedily merges adjacent triangle pairs into quads when the resulting quad
 * satisfies planarity, angle and shape quality criteria, then optionally
 * relaxes quad vertices toward planarity.
 */

#ifndef GEOMLIB_QUADREMESH_H
#define GEOMLIB_QUADREMESH_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float max_dihedral;        /* Max angle between the two triangle normals (deg) */
    float min_quad_angle;      /* Min interior quad angle (deg) */
    float max_quad_angle;      /* Max interior quad angle (deg) */
    int smoothing_iterations;  /* Planarity relaxation iterations */
} geomlib_quadremesh_params;

geomlib_quadremesh_params geomlib_quadremesh_params_default(void);

/* Mixed quad/triangle output mesh. Vertices are shared between quads and
 * remaining triangles. */
typedef struct {
    geomlib_vec3 *vertices;
    size_t vertex_count;
    uint32_t *quads;      /* 4 indices per quad */
    size_t quad_count;
    uint32_t *triangles;  /* 3 indices per remaining triangle */
    size_t triangle_count;
} geomlib_quadmesh;

void geomlib_quadmesh_free(geomlib_quadmesh *m);

/*
 * Convert a triangle mesh to a mixed quad/triangle mesh. Adjacent triangle
 * pairs that pass the quality criteria become quads; the rest stay triangles.
 */
int geomlib_remesh_to_quads(const geomlib_trimesh *input,
                            const geomlib_quadremesh_params *params,
                            geomlib_quadmesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_QUADREMESH_H */
