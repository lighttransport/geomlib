/*
 * geomlib - 2D segment intersection
 *
 * Robust classification and intersection of two planar segments.
 */

#ifndef GEOMLIB_SEGMENT2D_H
#define GEOMLIB_SEGMENT2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_SEG2D_NONE = 0,
    GEOMLIB_SEG2D_PROPER,      /* single interior crossing */
    GEOMLIB_SEG2D_ENDPOINT,    /* touch at an endpoint */
    GEOMLIB_SEG2D_COLLINEAR    /* overlap along a segment */
} geomlib_seg2d_kind;

typedef struct {
    geomlib_seg2d_kind kind;
    geomlib_vec2 point;        /* crossing point for PROPER/ENDPOINT */
    geomlib_vec2 overlap[2];   /* endpoints for COLLINEAR overlap */
    double s, t;               /* parameters on segment A and B */
} geomlib_seg2d_result;

geomlib_seg2d_result geomlib_segment2d_intersect(geomlib_vec2 a0, geomlib_vec2 a1,
                                                  geomlib_vec2 b0, geomlib_vec2 b1);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SEGMENT2D_H */
