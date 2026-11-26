/*
 * geomlib - 2D Delaunay refinement for minimum triangle angle
 *
 * Repeatedly inserts the circumcentre of the worst triangle when it lies
 * inside the domain, rebuilding a constrained Delaunay triangulation.
 */

#ifndef GEOMLIB_REFINE2D_H
#define GEOMLIB_REFINE2D_H

#include "geomlib_cdt2d.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float min_angle_deg;    /* Target minimum angle (default 25 degrees) */
    int max_points;         /* Steiner point budget (default 2000) */
} geomlib_refine2d_params;

geomlib_refine2d_params geomlib_refine2d_params_default(void);

/* Refine a simple polygon. Output vertices are owned by `out`. */
int geomlib_refine_polygon2d(const geomlib_vec2 *polygon, size_t n,
                             const geomlib_refine2d_params *params,
                             geomlib_cdt2d *out);

/* Smallest triangle angle in degrees. */
float geomlib_cdt2d_min_angle(const geomlib_cdt2d *c);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REFINE2D_H */
