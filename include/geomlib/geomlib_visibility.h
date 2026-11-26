/*
 * geomlib - Visibility graph shortest path among polygon obstacles
 *
 * Connects mutually visible obstacle vertices (and start/goal), then runs
 * Dijkstra to find the shortest collision-free path.
 */

#ifndef GEOMLIB_VISIBILITY_H
#define GEOMLIB_VISIBILITY_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *path;   /* start ... goal */
    size_t count;
    float length;
} geomlib_visibility_path;

geomlib_visibility_path geomlib_visibility_path_create(void);
void geomlib_visibility_path_free(geomlib_visibility_path *p);

/* `obstacles` is nh simple polygons with sizes `sizes`. All are CCW or CW. */
int geomlib_visibility_shortest_path(const geomlib_vec2 *start,
                                     const geomlib_vec2 *goal,
                                     const geomlib_vec2 *const *obstacles,
                                     const size_t *sizes, size_t nh,
                                     geomlib_visibility_path *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_VISIBILITY_H */
