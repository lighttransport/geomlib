/*
 * geomlib - Balanced mesh partition
 *
 * Assigns vertices to parts with balanced sizes using multi-source BFS region
 * growing from farthest-point seeds over the mesh adjacency graph.
 */

#ifndef GEOMLIB_PARTITION_H
#define GEOMLIB_PARTITION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *part_of_vertex;
    size_t vertex_count;
    size_t part_count;
    size_t *part_sizes;     /* Number of vertices per part */
} geomlib_partition_result;

geomlib_partition_result geomlib_partition_result_create(void);
void geomlib_partition_result_free(geomlib_partition_result *r);

/* Partition a mesh into `part_count` balanced regions. */
int geomlib_partition_mesh(const geomlib_trimesh *mesh,
                           size_t part_count,
                           geomlib_partition_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PARTITION_H */
