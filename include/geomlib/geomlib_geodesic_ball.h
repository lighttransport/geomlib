/*
 * geomlib - Geodesic ball extraction on triangle meshes
 *
 * Marks and lists all vertices within a geodesic (edge-graph) radius of a
 * source vertex.
 */

#ifndef GEOMLIB_GEODESIC_BALL_H
#define GEOMLIB_GEODESIC_BALL_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint8_t *mask;       /* one byte per vertex */
    uint32_t *vertices;  /* indices of included vertices */
    size_t count;
} geomlib_geodesic_ball;

geomlib_geodesic_ball geomlib_geodesic_ball_create(void);
void geomlib_geodesic_ball_free(geomlib_geodesic_ball *b);

int geomlib_extract_geodesic_ball(const geomlib_trimesh *mesh,
                                  uint32_t source, float radius,
                                  geomlib_geodesic_ball *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_BALL_H */
