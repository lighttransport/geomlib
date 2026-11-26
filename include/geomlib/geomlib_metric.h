/*
 * geomlib - Anisotropic metric tensors
 *
 * Symmetric positive-definite 3x3 metrics for size- and direction-dependent
 * distance control, with construction from principal curvatures. Intended for
 * metric-dependent discrete Voronoi clustering such as anisotropic ACVD.
 */

#ifndef GEOMLIB_METRIC_H
#define GEOMLIB_METRIC_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double m[3][3];   /* symmetric positive-definite matrix */
} geomlib_metric3;

geomlib_metric3 geomlib_metric_identity(void);
geomlib_metric3 geomlib_metric_diagonal(double l0, double l1, double l2);
geomlib_metric3 geomlib_metric_from_frame(const geomlib_vec3 axes[3],
                                          const double eig[3]);

double geomlib_metric_norm(geomlib_metric3 M, geomlib_vec3 v);
double geomlib_metric_distance(geomlib_metric3 M, geomlib_vec3 p,
                               geomlib_vec3 q);
double geomlib_metric_det(geomlib_metric3 M);
double geomlib_metric_anisotropy(geomlib_metric3 M);

/* Affine interpolation of two SPD metrics (the SPD cone is convex). */
geomlib_metric3 geomlib_metric_interpolate(geomlib_metric3 A,
                                           geomlib_metric3 B, double t);

/* Build a surface metric from principal curvatures. Tangent sizes shrink in
 * high-curvature directions; the normal size is scaled by `normal_scale`. */
int geomlib_metric_from_curvature(geomlib_vec3 normal, float k1, float k2,
                                  geomlib_vec3 d1, geomlib_vec3 d2,
                                  double target_size, double normal_scale,
                                  geomlib_metric3 *out);

/* Build one metric per mesh vertex from the mesh curvature field. `out` must
 * hold mesh->vertex_count entries. */
int geomlib_metric_field_from_curvature(const geomlib_trimesh *mesh,
                                        double target_size,
                                        double normal_scale,
                                        geomlib_metric3 *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_METRIC_H */
