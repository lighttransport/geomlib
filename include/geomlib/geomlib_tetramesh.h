/*
 * geomlib - Tetrahedral Mesh Generation for CFD
 * Delaunay tetrahedralization and volume meshing
 */

#ifndef GEOMLIB_TETRAMESH_H
#define GEOMLIB_TETRAMESH_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Tetrahedral Mesh Data Structure
 * ============================================================================ */

/* Tetrahedral element (4 vertices) */
typedef struct geomlib_tetrahedron {
    uint32_t vertices[4];  /* Indices to vertices */
} geomlib_tetrahedron;

/* Tetrahedral mesh structure */
typedef struct geomlib_tetramesh {
    geomlib_vec3 *vertices;           /* Vertex positions */
    geomlib_tetrahedron *tetrahedra;  /* Tetrahedral elements */

    size_t vertex_count;
    size_t tetrahedron_count;
    size_t vertex_capacity;
    size_t tetrahedron_capacity;

    /* Optional boundary information */
    uint32_t *boundary_faces;         /* Triangle indices on boundary (3 per face) */
    size_t boundary_face_count;
    bool *is_boundary_vertex;         /* Flags for boundary vertices */
} geomlib_tetramesh;

/* Tetrahedralization method */
typedef enum geomlib_tetra_method {
    GEOMLIB_TETRA_DELAUNAY,           /* Delaunay tetrahedralization */
    GEOMLIB_TETRA_ADVANCING_FRONT,    /* Advancing front method */
    GEOMLIB_TETRA_OCTREE              /* Octree-based meshing */
} geomlib_tetra_method;

/* Tetrahedralization parameters */
typedef struct geomlib_tetra_params {
    geomlib_tetra_method method;

    /* Size control */
    float max_edge_length;            /* Maximum edge length (0 = unlimited) */
    float min_edge_length;            /* Minimum edge length */
    float element_size;               /* Target element size */

    /* Quality control */
    float min_dihedral_angle;         /* Minimum dihedral angle (degrees) */
    float max_aspect_ratio;           /* Maximum aspect ratio */
    bool optimize_quality;            /* Apply mesh optimization */

    /* Boundary conforming */
    bool preserve_boundary;           /* Ensure boundary conformance */
    bool refine_boundary;             /* Refine near boundaries */

    /* Algorithm parameters */
    int max_iterations;               /* Max optimization iterations */
    int max_insertions;               /* Max Steiner point insertions */
} geomlib_tetra_params;

/* Tetrahedral mesh quality metrics */
typedef struct geomlib_tetra_quality {
    float volume;                     /* Element volume */
    float aspect_ratio;               /* Radius ratio (circumradius/inradius) */
    float min_dihedral_angle;         /* Minimum dihedral angle (degrees) */
    float max_dihedral_angle;         /* Maximum dihedral angle (degrees) */
    float edge_ratio;                 /* Longest/shortest edge ratio */
    float shape_measure;              /* Quality metric [0,1], 1=regular tetrahedron */
    bool is_sliver;                   /* Nearly degenerate element */
} geomlib_tetra_quality;

/* Mesh statistics */
typedef struct geomlib_tetramesh_stats {
    size_t vertex_count;
    size_t tetrahedron_count;
    size_t boundary_face_count;

    float total_volume;
    float min_volume;
    float max_volume;
    float avg_volume;

    float min_quality;
    float avg_quality;
    float worst_aspect_ratio;

    size_t sliver_count;
    size_t boundary_vertices;
} geomlib_tetramesh_stats;

/* ============================================================================
 * Tetrahedral Mesh Functions
 * ============================================================================ */

/* Initialize tetrahedral mesh */
void geomlib_tetramesh_init(
    geomlib_tetramesh *mesh,
    size_t vertex_capacity,
    size_t tetrahedron_capacity
);

/* Free tetrahedral mesh */
void geomlib_tetramesh_free(geomlib_tetramesh *mesh);

/* Add a vertex, returning its index or UINT32_MAX on allocation failure. */
uint32_t geomlib_tetramesh_add_vertex(
    geomlib_tetramesh *mesh,
    geomlib_vec3 position
);

/* Add a cell, returning its index or UINT32_MAX on allocation failure. */
uint32_t geomlib_tetramesh_add_tetrahedron(
    geomlib_tetramesh *mesh,
    uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3
);

/* Default parameters */
geomlib_tetra_params geomlib_tetra_params_default(geomlib_tetra_method method);

/* ============================================================================
 * Tetrahedralization Algorithms
 * ============================================================================ */

/* Fill a closed, consistently wound surface. A tetrahedral fan is used when
 * the volume centroid sees every face. Otherwise, a Delaunay candidate must
 * reproduce the boundary surface. Missing facets are recovered by inserting
 * Steiner points, so the boundary may be a subdivision of the input triangles;
 * the surface itself is preserved. Octree, boundary subdivision, and hard
 * maximum edge length are unsupported here. */
int geomlib_tetrahedralize_surface(
    const geomlib_trimesh *surface,
    const geomlib_tetra_params *params,
    geomlib_tetramesh *tetramesh
);

/* Delaunay tetrahedralization of a point cloud. The initialized output is
 * replaced on success. Coincident points retain their input vertex indices. */
int geomlib_tetrahedralize_points(
    const geomlib_vec3 *points,
    size_t point_count,
    geomlib_tetramesh *tetramesh
);

/* Generate tetrahedral mesh from bounding box */
int geomlib_tetrahedralize_box(
    geomlib_vec3 min_corner,
    geomlib_vec3 max_corner,
    const geomlib_tetra_params *params,
    geomlib_tetramesh *tetramesh
);

/* Preserve the boundary surface. A supplied interior point can serve as the fan
 * seed. The other interior points subdivide existing cells or are included in the
 * Delaunay candidate when the fan is unavailable. Missing constrained facets are
 * recovered by inserting Steiner points, so the boundary may be a subdivision of
 * the input triangles. The initialized output is replaced only on success. */
int geomlib_tetrahedralize_constrained(
    const geomlib_trimesh *boundary,
    const geomlib_vec3 *interior_points,
    size_t interior_count,
    geomlib_tetramesh *tetramesh
);

/* ============================================================================
 * Quality and Analysis
 * ============================================================================ */

/* Compute quality metrics for single tetrahedron */
geomlib_tetra_quality geomlib_compute_tetra_quality(
    const geomlib_tetramesh *mesh,
    uint32_t tetra_idx
);

/* Compute mesh statistics */
void geomlib_tetramesh_compute_stats(
    const geomlib_tetramesh *mesh,
    geomlib_tetramesh_stats *stats
);

/* Compute tetrahedron volume */
float geomlib_tetrahedron_volume(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
);

/* Interior dihedral angle in degrees for edges 01, 02, 03, 12, 13, 23.
 * Returns NAN when the edge or adjacent face is degenerate. */
float geomlib_tetrahedron_dihedral_angle(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3,
    int edge_idx  /* Edge index 0-5 */
);

/* Check if point is inside tetrahedron */
bool geomlib_point_in_tetrahedron(
    geomlib_vec3 point,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3
);

/* Compute circumsphere of tetrahedron */
void geomlib_tetrahedron_circumsphere(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_vec3 v3,
    geomlib_vec3 *center,
    float *radius
);

/* ============================================================================
 * Mesh Optimization
 * ============================================================================ */

/* Attempt up to iterations local interior flips (2-to-3 face flips and 3-to-2
 * edge flips) to improve worst-cell shape. Boundary triangles stay fixed. */
int geomlib_tetramesh_optimize(
    geomlib_tetramesh *mesh,
    int iterations
);

/* Remove slivers where a valid interior-face flip reduces the local sliver count.
 * A successful pass can retain slivers that require other operations. */
int geomlib_tetramesh_remove_slivers(
    geomlib_tetramesh *mesh
);

/* Laplacian smoothing of interior vertices. Boundary vertices and boundary
 * triangles stay fixed. A move is accepted only when every incident cell keeps
 * positive volume and local worst-cell quality improves. */
int geomlib_tetramesh_smooth(
    geomlib_tetramesh *mesh,
    int iterations
);

/* Bisect the longest eligible edge across all incident cells. The target is
 * max_edge_length, or element_size when max_edge_length is zero. The chosen
 * edge must exceed twice min_edge_length. Both budgets limit edge bisections.
 * Boundary edges are locked unless refine_boundary is true and preserve_boundary
 * is false. A successful pass returns the mesh reached within those constraints. */
int geomlib_tetramesh_refine(
    geomlib_tetramesh *mesh,
    const geomlib_tetra_params *params
);

/* ============================================================================
 * Boundary Operations
 * ============================================================================ */

/* Append outward boundary triangles and the source vertices to an initialized mesh. */
int geomlib_tetramesh_extract_boundary(
    const geomlib_tetramesh *tetramesh,
    geomlib_trimesh *boundary
);

/* Refresh outward boundary faces and boundary vertex flags on success. */
void geomlib_tetramesh_mark_boundary(
    geomlib_tetramesh *mesh
);

/* ============================================================================
 * Export Functions
 * ============================================================================ */

/* Export an ASCII VTK legacy unstructured grid or an ASCII Gmsh 2.2 mesh. */
int geomlib_tetramesh_export_vtk(
    const geomlib_tetramesh *mesh,
    const char *filename
);

int geomlib_tetramesh_export_gmsh(
    const geomlib_tetramesh *mesh,
    const char *filename
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TETRAMESH_H */
