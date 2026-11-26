/*
 * geomlib - 2D polygon Boolean operations
 *
 * Splits both polygon boundaries at their crossings, classifies each sub-edge
 * by containment in the other polygon, keeps the relevant parts, and traces
 * the result into closed loops.
 */

#ifndef GEOMLIB_POLYBOOL2D_H
#define GEOMLIB_POLYBOOL2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_POLY_UNION = 0,
    GEOMLIB_POLY_INTERSECTION = 1,
    GEOMLIB_POLY_DIFFERENCE = 2   /* A minus B */
} geomlib_poly_bool_op;

typedef struct {
    geomlib_vec2 *loops;     /* concatenated loop vertices */
    size_t *offsets;         /* loop_count+1 offsets into loops */
    size_t loop_count;
} geomlib_poly_bool_result;

geomlib_poly_bool_result geomlib_poly_bool_result_create(void);
void geomlib_poly_bool_result_free(geomlib_poly_bool_result *r);

/* Both inputs are simple polygons in any winding. Returns GEOMLIB_MESH_SUCCESS. */
int geomlib_polygon_boolean(const geomlib_vec2 *a, size_t na,
                            const geomlib_vec2 *b, size_t nb,
                            geomlib_poly_bool_op op,
                            geomlib_poly_bool_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYBOOL2D_H */
