/*
 * geomlib - Mesh Smoothing
 * Laplacian and Taubin smoothing algorithms
 */

#ifndef GEOMLIB_SMOOTHING_H
#define GEOMLIB_SMOOTHING_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Smoothing method */
typedef enum geomlib_smoothing_method {
    GEOMLIB_SMOOTH_LAPLACIAN,    /* Basic Laplacian smoothing (may cause shrinkage) */
    GEOMLIB_SMOOTH_TAUBIN,       /* Taubin lambda/mu pair per iteration (volume preserving) */
    GEOMLIB_SMOOTH_HC            /* HC (Humphrey's Classes) feature correction */
} geomlib_smoothing_method;

/* Smoothing parameters */
typedef struct geomlib_smoothing_params {
    geomlib_smoothing_method method;

    /* Laplacian parameters */
    float lambda;                /* Smoothing weight (0.0 to 1.0, typical: 0.5) */
    int iterations;              /* Number of smoothing iterations */

    /* Taubin parameters */
    float mu;                    /* Second smoothing weight (negative, typical: -0.53) */

    /* HC parameters */
    float alpha;                 /* HC alpha parameter (typical: 0.5) */
    float beta;                  /* HC beta parameter (typical: 0.5) */

    /* General */
    bool preserve_boundaries;    /* Don't move boundary vertices */
    bool preserve_features;      /* Freeze vertices incident to detected sharp edges */
    float feature_angle;         /* Sharp-edge threshold in degrees */
} geomlib_smoothing_params;

/* Create default smoothing parameters */
geomlib_smoothing_params geomlib_smoothing_params_default(geomlib_smoothing_method method);

/* Smooth mesh using specified method */
int geomlib_smooth_mesh(
    geomlib_trimesh *mesh,
    const geomlib_smoothing_params *params
);

/* Simplified interfaces */

/* Laplacian smoothing (simple, may cause shrinkage) */
int geomlib_smooth_laplacian(
    geomlib_trimesh *mesh,
    float lambda,
    int iterations
);

/* Taubin smoothing (preserves volume better) */
int geomlib_smooth_taubin(
    geomlib_trimesh *mesh,
    float lambda,
    float mu,
    int iterations
);

/* Smooth only vertex normals (for shading) */
int geomlib_smooth_normals(
    geomlib_trimesh *mesh,
    int iterations
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SMOOTHING_H */
