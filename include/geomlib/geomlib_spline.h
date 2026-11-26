/*
 * geomlib - Catmull-Rom spline resampling
 *
 * Resamples a 3D polyline by evaluating uniform or centripetal Catmull-Rom
 * segments with the Barry-Goldman pyramid.
 */

#ifndef GEOMLIB_SPLINE_H
#define GEOMLIB_SPLINE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_CATMULL_UNIFORM = 0,
    GEOMLIB_CATMULL_CENTRIPETAL = 1
} geomlib_catmull_type;

/* Resample a polyline into `out`. For open curves the output length is
 * (n-1)*samples_per_segment + 1; for closed curves it is n*samples_per_segment
 * plus a repeated first point (total n*samples_per_segment + 1). */
int geomlib_resample_catmull_rom(const geomlib_vec3 *points, size_t n,
                                 int closed, size_t samples_per_segment,
                                 geomlib_catmull_type type,
                                 geomlib_vec3 *out);

/* Number of output points produced by geomlib_resample_catmull_rom. */
size_t geomlib_catmull_output_count(size_t n, int closed,
                                    size_t samples_per_segment);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPLINE_H */
