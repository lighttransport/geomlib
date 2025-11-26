/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Monte Carlo Geometry Processing
 * Core functions implemented in C for maximum portability
 *
 * References:
 * - Keenan Crane et al. "Monte Carlo Geometry Processing" (SIGGRAPH 2020)
 * - Matt Pharr "Adventures in Sampling Points on Triangles"
 * - "Physically Based Rendering" pbr-book.org
 */

#ifndef GEOMLIB_MONTECARLO_H
#define GEOMLIB_MONTECARLO_H

#include <stddef.h>
#include <stdint.h>
#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include "geomlib_random.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Triangle Surface Sampling
 * ============================================================================ */

/**
 * Sample a uniform random point on a triangle using barycentric coordinates
 *
 * Uses the standard transformation: b0 = 1 - sqrt(u), b1 = v * sqrt(u)
 * This ensures uniform distribution across the triangle area.
 *
 * @param v0        First vertex of triangle
 * @param v1        Second vertex of triangle
 * @param v2        Third vertex of triangle
 * @param rng       PCG32 random number generator
 * @return          Random point on triangle surface
 */
geomlib_vec3 geomlib_sample_triangle_uniform(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_pcg32 *rng
);

/**
 * Sample a uniform random point on a triangle with barycentric coordinates
 *
 * Same as geomlib_sample_triangle_uniform but also returns the barycentric
 * coordinates used, which can be useful for interpolating attributes.
 *
 * @param v0        First vertex of triangle
 * @param v1        Second vertex of triangle
 * @param v2        Third vertex of triangle
 * @param rng       PCG32 random number generator
 * @param bary_out  Output barycentric coordinates (b0, b1, b2) where b0+b1+b2=1
 * @return          Random point on triangle surface
 */
geomlib_vec3 geomlib_sample_triangle_barycentric(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_pcg32 *rng,
    geomlib_vec3 *bary_out
);

/**
 * Compute the PDF for uniform triangle sampling
 *
 * For uniform sampling, PDF = 1 / area
 *
 * @param v0        First vertex of triangle
 * @param v1        Second vertex of triangle
 * @param v2        Third vertex of triangle
 * @return          PDF value (1 / triangle_area)
 */
float geomlib_triangle_sample_pdf(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
);

/* ============================================================================
 * Mesh Surface Sampling
 * ============================================================================ */

/**
 * Sample distribution for mesh surface sampling
 *
 * Stores cumulative distribution function (CDF) for area-weighted sampling.
 * This allows efficient O(log n) sampling from the mesh surface.
 */
typedef struct geomlib_mesh_sample_distribution {
    float *triangle_areas;      /* Area of each triangle */
    float *cdf;                 /* Cumulative distribution function */
    float total_area;           /* Total surface area */
    size_t triangle_count;      /* Number of triangles */
} geomlib_mesh_sample_distribution;

/**
 * Build sampling distribution for a triangle mesh
 *
 * Computes the CDF for area-weighted sampling. This is a preprocessing
 * step that enables fast O(log n) sampling.
 *
 * @param mesh      Triangle mesh to build distribution for
 * @param dist      Distribution structure to initialize
 * @return          GEOMLIB_MESH_SUCCESS on success, error code otherwise
 */
int geomlib_mesh_build_sample_distribution(
    const geomlib_trimesh *mesh,
    geomlib_mesh_sample_distribution *dist
);

/**
 * Free sampling distribution resources
 *
 * @param dist      Distribution to free
 */
void geomlib_mesh_free_sample_distribution(
    geomlib_mesh_sample_distribution *dist
);

/**
 * Sample a uniform random point on mesh surface
 *
 * Uses area-weighted sampling to ensure uniform distribution across
 * the entire mesh surface. Requires a precomputed distribution.
 *
 * @param mesh      Triangle mesh to sample from
 * @param dist      Precomputed sampling distribution
 * @param rng       PCG32 random number generator
 * @param triangle_idx_out  Optional output for selected triangle index (can be NULL)
 * @return          Random point on mesh surface
 */
geomlib_vec3 geomlib_sample_mesh_surface(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_sample_distribution *dist,
    geomlib_pcg32 *rng,
    size_t *triangle_idx_out
);

/**
 * Sample a uniform random point on mesh surface with normal
 *
 * Same as geomlib_sample_mesh_surface but also returns the interpolated
 * normal at the sample point.
 *
 * @param mesh      Triangle mesh to sample from (must have normals computed)
 * @param dist      Precomputed sampling distribution
 * @param rng       PCG32 random number generator
 * @param normal_out Output interpolated normal at sample point
 * @param triangle_idx_out  Optional output for selected triangle index (can be NULL)
 * @return          Random point on mesh surface
 */
geomlib_vec3 geomlib_sample_mesh_surface_with_normal(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_sample_distribution *dist,
    geomlib_pcg32 *rng,
    geomlib_vec3 *normal_out,
    size_t *triangle_idx_out
);

/* ============================================================================
 * Monte Carlo Integration Helpers
 * ============================================================================ */

/**
 * Monte Carlo estimator for surface integral
 *
 * Estimates the integral: ∫_S f(x) dA
 * using N samples: (1/N) * Σ f(x_i) / pdf(x_i)
 *
 * For uniform sampling on surface, pdf(x) = 1/area, so:
 * Estimate = (area/N) * Σ f(x_i)
 *
 * @param samples       Array of sample points
 * @param values        Function values f(x_i) at each sample
 * @param pdf           PDF values at each sample (or NULL for uniform with area)
 * @param num_samples   Number of samples
 * @param area          Total surface area (used if pdf is NULL)
 * @return              Estimated integral value
 */
double geomlib_mc_estimate_surface_integral(
    const geomlib_vec3 *samples,
    const float *values,
    const float *pdf,
    size_t num_samples,
    float area
);

/**
 * Compute variance of Monte Carlo estimator
 *
 * Estimates the variance: Var[estimator] = (1/N) * Var[f(X)/pdf(X)]
 *
 * @param values        Function values f(x_i) at each sample
 * @param pdf           PDF values at each sample
 * @param num_samples   Number of samples
 * @param mean          Mean value (can be precomputed or pass estimate)
 * @return              Estimated variance
 */
double geomlib_mc_estimate_variance(
    const float *values,
    const float *pdf,
    size_t num_samples,
    double mean
);

/* ============================================================================
 * Closest Point Queries (for Walk on Spheres)
 * ============================================================================ */

/**
 * Find closest point on triangle to query point
 *
 * Returns the closest point on the triangle surface (including edges and vertices)
 * to the given query point.
 *
 * @param query     Query point
 * @param v0        First vertex of triangle
 * @param v1        Second vertex of triangle
 * @param v2        Third vertex of triangle
 * @return          Closest point on triangle to query
 */
geomlib_vec3 geomlib_closest_point_triangle(
    geomlib_vec3 query,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
);

/**
 * Find closest point on mesh to query point
 *
 * Brute force search through all triangles. For production use,
 * consider building a spatial acceleration structure.
 *
 * @param query         Query point
 * @param mesh          Triangle mesh
 * @param distance_out  Optional output for distance to closest point (can be NULL)
 * @return              Closest point on mesh to query
 */
geomlib_vec3 geomlib_closest_point_mesh(
    geomlib_vec3 query,
    const geomlib_trimesh *mesh,
    float *distance_out
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Sample uniform random point on unit sphere surface
 *
 * Uses the standard spherical coordinate method.
 *
 * @param rng       PCG32 random number generator
 * @return          Random point on unit sphere
 */
geomlib_vec3 geomlib_sample_sphere_uniform(geomlib_pcg32 *rng);

/**
 * Sample uniform random point in unit sphere volume
 *
 * Uses rejection sampling in the unit cube.
 *
 * @param rng       PCG32 random number generator
 * @return          Random point inside unit sphere
 */
geomlib_vec3 geomlib_sample_sphere_volume(geomlib_pcg32 *rng);

/**
 * Sample uniform random direction in hemisphere
 *
 * Generates random direction in hemisphere aligned with given normal.
 * Useful for Monte Carlo rendering and geometry processing.
 *
 * @param normal    Normal vector defining hemisphere orientation
 * @param rng       PCG32 random number generator
 * @return          Random direction in hemisphere
 */
geomlib_vec3 geomlib_sample_hemisphere_uniform(geomlib_vec3 normal, geomlib_pcg32 *rng);

/**
 * Sample cosine-weighted direction in hemisphere
 *
 * Generates random direction with pdf proportional to cos(theta).
 * Commonly used in diffuse reflection Monte Carlo integration.
 *
 * @param normal    Normal vector defining hemisphere orientation
 * @param rng       PCG32 random number generator
 * @return          Random direction in hemisphere (cosine weighted)
 */
geomlib_vec3 geomlib_sample_hemisphere_cosine(geomlib_vec3 normal, geomlib_pcg32 *rng);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MONTECARLO_H */
