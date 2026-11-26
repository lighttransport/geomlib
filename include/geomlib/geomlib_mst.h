/*
 * geomlib - Minimum spanning tree of a point set
 *
 * Prim's algorithm on the complete Euclidean graph.
 */

#ifndef GEOMLIB_MST_H
#define GEOMLIB_MST_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *edges;   /* 2 indices per MST edge */
    size_t edge_count;
    double total_weight;
} geomlib_mst_result;

geomlib_mst_result geomlib_mst_result_create(void);
void geomlib_mst_result_free(geomlib_mst_result *r);

/* Compute the MST of `count` points (dim coordinates per point, 2 or 3). */
int geomlib_compute_mst(const double *points, size_t count, int dim,
                        geomlib_mst_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MST_H */
