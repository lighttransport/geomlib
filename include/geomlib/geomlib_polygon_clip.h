/*
 * geomlib - Sutherland-Hodgman convex polygon clipping
 *
 * Clips a subject polygon against a convex clipping polygon.
 */

#ifndef GEOMLIB_POLYGON_CLIP_H
#define GEOMLIB_POLYGON_CLIP_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Clip `subject` by convex `clip` (both CCW) into the caller's `geomlib_polygon2d`. */
int geomlib_clip_polygon_convex(const geomlib_vec2 *subject, size_t ns,
                                const geomlib_vec2 *clip, size_t nc,
                                geomlib_polygon2d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_CLIP_H */
