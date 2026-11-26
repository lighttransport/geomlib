/*
 * geomlib - Broad-phase triangle overlap queries
 *
 * Builds bounding-box trees over two triangle meshes and reports all pairs of
 * triangles whose AABBs overlap.
 */

#ifndef GEOMLIB_MESH_OVERLAP_H
#define GEOMLIB_MESH_OVERLAP_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *pairs;    /* two triangle indices per overlapping pair */
    size_t count;
    size_t capacity;
} geomlib_overlap_pairs;

geomlib_overlap_pairs geomlib_overlap_pairs_create(void);
void geomlib_overlap_pairs_free(geomlib_overlap_pairs *p);

/* Report all AABB-overlapping triangle pairs between meshes a and b. */
int geomlib_aabb_overlap_pairs(const geomlib_trimesh *a,
                               const geomlib_trimesh *b,
                               geomlib_overlap_pairs *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MESH_OVERLAP_H */
