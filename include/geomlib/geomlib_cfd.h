/*
 * geomlib - CFD (Computational Fluid Dynamics) Utilities
 * Gradient computation, interpolation, streamlines, and mesh quality for CFD
 */

#ifndef GEOMLIB_CFD_H
#define GEOMLIB_CFD_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include "geomlib_tetramesh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Scalar and Vector Fields
 * ============================================================================ */

/* Scalar field (e.g., pressure, temperature) */
typedef struct geomlib_scalar_field {
    float *values;              /* Field values at vertices/cells */
    size_t count;
    geomlib_vec3 *gradients;    /* Optional: gradient at each point */
} geomlib_scalar_field;

/* Vector field (e.g., velocity, force) */
typedef struct geomlib_vector_field {
    geomlib_vec3 *values;       /* Field values at vertices/cells */
    size_t count;
} geomlib_vector_field;

/* Initialize fields */
geomlib_scalar_field geomlib_scalar_field_create(size_t count);
void geomlib_scalar_field_free(geomlib_scalar_field *field);
geomlib_vector_field geomlib_vector_field_create(size_t count);
void geomlib_vector_field_free(geomlib_vector_field *field);

/* ============================================================================
 * Face-Based Mesh for Finite Volume Method
 * ============================================================================ */

/* Face structure for FVM */
typedef struct geomlib_fvm_face {
    uint32_t vertices[4];       /* Vertex indices (up to quad, use 3 for triangle) */
    uint32_t vertex_count;      /* 3 for triangle, 4 for quad */
    uint32_t cell_left;         /* Left cell index */
    uint32_t cell_right;        /* Right cell index (UINT32_MAX for boundary) */
    geomlib_vec3 center;        /* Face center */
    geomlib_vec3 normal;        /* Face normal (pointing left->right) */
    float area;                 /* Face area */
} geomlib_fvm_face;

/* Face-based mesh structure */
typedef struct geomlib_fvm_mesh {
    /* Original mesh */
    geomlib_tetramesh *tetramesh;  /* Borrowed source volume mesh */

    /* Face data */
    geomlib_fvm_face *faces;
    size_t face_count;

    /* Cell data */
    geomlib_vec3 *cell_centers;    /* Cell centroids */
    float *cell_volumes;           /* Cell volumes */
    size_t cell_count;

    /* Connectivity */
    uint32_t **cell_faces;         /* Faces adjacent to each cell */
    uint32_t *cell_face_count;     /* Number of faces per cell */
} geomlib_fvm_mesh;

/* Build unique faces and cell geometry. The source must outlive the result.
 * Use an empty output; it is assigned only on success. */
int geomlib_fvm_build(
    geomlib_fvm_mesh *fvm,
    const geomlib_tetramesh *tetramesh
);

/* Free FVM mesh */
void geomlib_fvm_free(geomlib_fvm_mesh *fvm);

/* ============================================================================
 * Gradient Computation
 * ============================================================================ */

/* Gradient method */
typedef enum geomlib_gradient_method {
    GEOMLIB_GRADIENT_GREEN_GAUSS,      /* Green-Gauss cell-based */
    GEOMLIB_GRADIENT_GREEN_GAUSS_NODE, /* Green-Gauss node-based */
    GEOMLIB_GRADIENT_LEAST_SQUARES     /* Least-squares reconstruction */
} geomlib_gradient_method;

/* Allocate the vertex gradients of a nodal scalar field on a triangle mesh.
 * Each method evaluates the same linear field on each triangle, then averages
 * the tangential gradients at incident vertices using triangle area. */
int geomlib_compute_gradient_scalar(
    const geomlib_trimesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_gradient_method method,
    geomlib_vector_field *gradient
);

/* Allocate one gradient per cell from scalar values at the tetrahedral vertices.
 * Each supported method recovers the gradient of the linear interpolant. */
int geomlib_compute_gradient_tetra(
    const geomlib_tetramesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_gradient_method method,
    geomlib_vector_field *gradient
);

/* Green-Gauss gradient of a planar polygon with linearly varying values along each edge. */
geomlib_vec3 geomlib_gradient_green_gauss(
    const geomlib_vec3 *vertices,
    const float *values,
    size_t count
);

/* Three-dimensional least-squares fit; neighbor_values are offsets from the center value. */
geomlib_vec3 geomlib_gradient_least_squares(
    geomlib_vec3 center,
    const geomlib_vec3 *neighbor_positions,
    const float *neighbor_values,
    size_t neighbor_count
);

/* ============================================================================
 * Field Operations
 * ============================================================================ */

/* Tangential trace of the Jacobian of an ambient nodal vector field, averaged by area. */
int geomlib_compute_divergence(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *field,
    geomlib_scalar_field *divergence
);

/* Normal component of the nodal vector field curl, averaged by incident triangle area. */
int geomlib_compute_curl(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *field,
    geomlib_vector_field *curl
);

/* Cotangent Laplace-Beltrami operator with lumped vertex areas and the weak
 * natural boundary convention. Field operators write to preallocated outputs. */
int geomlib_compute_laplacian(
    const geomlib_trimesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_scalar_field *laplacian
);

/* ============================================================================
 * Interpolation Schemes
 * ============================================================================ */

/* Interpolation scheme */
typedef enum geomlib_interpolation_scheme {
    GEOMLIB_INTERP_LINEAR,          /* Linear interpolation */
    GEOMLIB_INTERP_UPWIND,          /* Upwind (first-order) */
    GEOMLIB_INTERP_CENTRAL,         /* Central differencing */
    GEOMLIB_INTERP_QUICK,           /* QUICK scheme */
    GEOMLIB_INTERP_MUSCL            /* MUSCL scheme */
} geomlib_interpolation_scheme;

/* Sample the nodal field at the closest surface point using barycentric weights.
 * Point sampling uses linear weights for all scheme values. */
float geomlib_interpolate_scalar(
    const geomlib_trimesh *mesh,
    const geomlib_scalar_field *field,
    geomlib_vec3 point,
    geomlib_interpolation_scheme scheme
);

/* Vector version of geomlib_interpolate_scalar. */
geomlib_vec3 geomlib_interpolate_vector(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *field,
    geomlib_vec3 point,
    geomlib_interpolation_scheme scheme
);

/* Face value interpolation for FVM */
float geomlib_interpolate_face_value(
    float value_left,
    float value_right,
    float distance_ratio,
    geomlib_interpolation_scheme scheme
);

/* ============================================================================
 * Streamline Tracing
 * ============================================================================ */

/* Streamline integration method */
typedef enum geomlib_integration_method {
    GEOMLIB_INTEGRATE_EULER,        /* Forward Euler */
    GEOMLIB_INTEGRATE_RK2,          /* Runge-Kutta 2nd order */
    GEOMLIB_INTEGRATE_RK4           /* Runge-Kutta 4th order */
} geomlib_integration_method;

/* Streamline structure */
typedef struct geomlib_streamline {
    geomlib_vec3 *points;           /* Points along streamline */
    size_t point_count;
    size_t point_capacity;
} geomlib_streamline;

/* Streamline parameters */
typedef struct geomlib_streamline_params {
    geomlib_integration_method method;
    float step_size;                /* Integration step size */
    int max_steps;                  /* Maximum number of steps */
    float max_length;               /* Maximum streamline length */
    geomlib_aabb bounds;            /* Bounding box (stop if exited) */
} geomlib_streamline_params;

/* Initialize streamline */
geomlib_streamline geomlib_streamline_create(size_t capacity);

/* Free streamline */
void geomlib_streamline_free(geomlib_streamline *streamline);

/* Trace streamline from seed point */
int geomlib_trace_streamline(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    geomlib_vec3 seed,
    const geomlib_streamline_params *params,
    geomlib_streamline *streamline
);

/* Generate multiple streamlines */
int geomlib_generate_streamlines(
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    const geomlib_vec3 *seeds,
    size_t seed_count,
    const geomlib_streamline_params *params,
    geomlib_streamline *streamlines  /* Array of seed_count streamlines */
);

/* ============================================================================
 * CFD Mesh Quality
 * ============================================================================ */

/* CFD-specific quality metrics */
typedef struct geomlib_cfd_quality {
    float skewness;                 /* Edge-length coefficient of variation, clamped to [0,1] */
    float non_orthogonality;        /* Maximum interior-face non-orthogonality angle, in degrees */
    float expansion_ratio;          /* Maximum larger-to-smaller volume ratio across interior faces */
    float face_flatness;            /* Face warping metric */
    bool is_valid;                  /* Meets CFD quality criteria */
} geomlib_cfd_quality;

/* CFD quality thresholds */
typedef struct geomlib_cfd_quality_thresholds {
    float max_skewness;             /* Default: 0.95 */
    float max_non_orthogonality;    /* Default: 70 degrees */
    float max_expansion_ratio;      /* Default: 5.0 */
} geomlib_cfd_quality_thresholds;

/* Compute CFD quality for tetrahedron */
geomlib_cfd_quality geomlib_compute_cfd_quality(
    const geomlib_tetramesh *mesh,
    uint32_t tet_idx,
    const geomlib_cfd_quality_thresholds *thresholds
);

/* Check mesh quality for CFD */
bool geomlib_check_cfd_mesh_quality(
    const geomlib_tetramesh *mesh,
    const geomlib_cfd_quality_thresholds *thresholds,
    size_t *failed_count
);

/* Default quality thresholds */
geomlib_cfd_quality_thresholds geomlib_cfd_quality_thresholds_default(void);

/* ============================================================================
 * Particle Tracing
 * ============================================================================ */

/* Particle structure */
typedef struct geomlib_particle {
    geomlib_vec3 position;
    geomlib_vec3 velocity;
    float time;
    bool active;
} geomlib_particle;

/* Particle system */
typedef struct geomlib_particle_system {
    geomlib_particle *particles;
    size_t particle_count;
    size_t particle_capacity;
} geomlib_particle_system;

/* Create particle system */
geomlib_particle_system geomlib_particle_system_create(size_t capacity);

/* Free particle system */
void geomlib_particle_system_free(geomlib_particle_system *system);

/* Advect particles in velocity field */
int geomlib_advect_particles(
    geomlib_particle_system *particles,
    const geomlib_trimesh *mesh,
    const geomlib_vector_field *velocity,
    float dt,
    const geomlib_aabb *bounds
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/* Compute cell center (tetrahedron centroid) */
geomlib_vec3 geomlib_tetrahedron_centroid(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
);

/* Compute face area and normal for triangle */
void geomlib_face_area_normal(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float *area,
    geomlib_vec3 *normal
);

/* Distance from point to plane */
float geomlib_point_plane_distance(
    geomlib_vec3 point,
    geomlib_vec3 plane_point,
    geomlib_vec3 plane_normal
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CFD_H */
