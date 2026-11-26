/*
 * geomlib - Feature-preserving point-set processing
 *
 * Bilateral denoising (smoothing constrained by normal similarity) and
 * statistical outlier removal for unorganized point clouds.
 */

#ifndef GEOMLIB_POINTSET_H
#define GEOMLIB_POINTSET_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float spatial_sigma;   /* Spatial Gaussian width (default: average spacing) */
    float normal_sigma;    /* Normal-similarity width (default 0.5) */
    int neighbors;         /* k nearest neighbours (default 16) */
    int iterations;        /* Denoising iterations (default 2) */
} geomlib_bilateral_params;

geomlib_bilateral_params geomlib_bilateral_params_default(void);

/* Smooth point positions in place, preserving edges via normal weighting. */
int geomlib_pointcloud_bilateral_smooth(geomlib_pointcloud *cloud,
                                        const geomlib_bilateral_params *params);

/* Remove the fraction of points with the largest average k-NN distance.
 * Returns the number of removed points. */
size_t geomlib_pointcloud_remove_outliers(geomlib_pointcloud *cloud,
                                          float fraction, int k);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POINTSET_H */
