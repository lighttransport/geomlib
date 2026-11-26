# geomlib

A portable geometry processing library in C++17 with core functions implemented in C.

## Features

- **Linear Algebra**: Vectors (vec2, vec3, vec4), matrices (mat3, mat4), transformations
- **Linear Solvers**: Gaussian elimination, LU decomposition, Conjugate Gradient, Jacobi iteration
- **Mesh Processing**: Triangle meshes, polygon meshes, primitive generation
- **Geometry Processing**: Heat-method geodesics, HC/Taubin smoothing, Loop and modified Butterfly subdivision, point-cloud reconstruction
- **No External Dependencies**: Custom linear algebra implementation (no Eigen required)
- **C++17 Compatible**: No RTTI, no C++ exceptions
- **Cross-platform**: Core functions in C for maximum portability

## Building

```bash
mkdir build && cd build
cmake ..
make
```

### Build Options

- `GEOMLIB_BUILD_TESTS` - Build tests (default: ON)
- `GEOMLIB_BUILD_EXAMPLES` - Build examples (default: OFF)

## Usage

```cpp
#include <geomlib/geomlib.h>

// Linear algebra
geomlib_vec3 a = geomlib_vec3_create(1.0f, 2.0f, 3.0f);
geomlib_vec3 b = geomlib_vec3_create(4.0f, 5.0f, 6.0f);
geomlib_vec3 c = geomlib_vec3_cross(a, b);
float dot = geomlib_vec3_dot(a, b);

// Matrix operations
geomlib_mat4 transform = geomlib_mat4_translate(1.0f, 2.0f, 3.0f);
geomlib_mat4 rotation = geomlib_mat4_rotate_y(3.14159f / 4.0f);
geomlib_mat4 combined = geomlib_mat4_mul(transform, rotation);

// Mesh creation
geomlib_trimesh mesh;
geomlib_create_sphere(&mesh, 1.0f, 32, 16);
geomlib_trimesh_compute_normals(&mesh);
geomlib_aabb bounds = geomlib_trimesh_compute_aabb(&mesh);
geomlib_trimesh_free(&mesh);

// Linear solver
float A[] = {4.0f, 1.0f, 1.0f, 3.0f};
float b[] = {1.0f, 2.0f};
geomlib_solve_gauss(2, A, b);  // Solution in b
```

## API Overview

### Linear Algebra (`geomlib_linalg.h`)

- Vector types: `geomlib_vec2`, `geomlib_vec3`, `geomlib_vec4` (float and double variants)
- Matrix types: `geomlib_mat3`, `geomlib_mat4` (column-major order)
- Vector operations: add, sub, mul, div, dot, cross, length, normalize
- Matrix operations: multiply, transpose, determinant, inverse
- Transformation matrices: translate, scale, rotate, look_at, perspective, ortho

### Solvers (`geomlib_solver.h`)

- `geomlib_solve_gauss()` - Gaussian elimination with partial pivoting
- `geomlib_lu_decompose()` / `geomlib_lu_solve()` - LU decomposition
- `geomlib_solve_cg()` - Conjugate Gradient for symmetric positive-definite systems
- `geomlib_solve_jacobi()` - Jacobi iterative solver
- `geomlib_eigen_power()` - Power iteration for dominant eigenvalue
- `geomlib_solve_least_squares()` - Least squares solver

### Mesh (`geomlib_mesh.h`)

- Polygon types: `geomlib_polygon2d`, `geomlib_polygon3d`
- Mesh types: `geomlib_trimesh`, `geomlib_polymesh`
- Bounding box: `geomlib_aabb`
- Polygon operations: area, perimeter, centroid, convexity test, point containment
- Mesh operations: add vertices/triangles, compute normals, compute bounds
- Primitive generation: box, sphere, cylinder, plane
- Conversion: triangulate polygon mesh

## Topology, validation, and SDF

Use `geomlib_find_connected_components` for component labels and
`geomlib_extract_boundary_loops` for closed boundary cycles. Release the
results with `geomlib_component_info_free` and `geomlib_boundary_loops_free`.

Run Basic validation before processing mesh data from an external source.
Standard also checks vertex coordinates, unused vertices, manifoldness, and
winding. Strict adds self-intersection checks and rejects near-zero-area
triangles. An open boundary can be valid. Use `geomlib_analyze_manifold` to
check `is_closed` before generating a signed distance field. Use unsigned
generation for open surfaces by setting `compute_sign` to `false`.

For edge, angle, and area summaries, call `geomlib_compute_mesh_stats(&mesh, 16)`
with the desired number of histogram bins. The returned struct owns the
histograms. Release them with `geomlib_mesh_stats_free`.

```c
/* mesh is an existing triangle mesh. */
geomlib_validation_result check = geomlib_validate_mesh(&mesh, GEOMLIB_VALIDATE_STRICT);
geomlib_manifold_analysis topology = {0};
if (check.is_valid) topology = geomlib_analyze_manifold(&mesh);
if (topology.is_closed) {
    geomlib_sdf sdf = geomlib_mesh_to_sdf_resolution(&mesh, 32);
    geomlib_trimesh rebuilt = {0};
    if (sdf.values && geomlib_sdf_extract_surface(&sdf, &rebuilt) == GEOMLIB_MESH_SUCCESS) {
        /* Use rebuilt, then release its vertex and index arrays. */
        geomlib_trimesh_free(&rebuilt);
    }
    geomlib_sdf_free(&sdf);
}
geomlib_manifold_analysis_free(&topology);
geomlib_validation_result_free(&check);
```

SDF queries and extraction use voxel-center samples. The SDF and reconstructed
mesh own separate storage, so either can be freed first. See the
[Phase 1 and 2 roadmap](PHASE_1_2_ROADMAP.md) for current coverage and follow-up work.

Generation defaults to `GEOMLIB_SDF_ALGORITHM_BRUTE_FORCE`, which evaluates every
voxel against a BVH. Set `algorithm = GEOMLIB_SDF_ALGORITHM_FAST_SWEEPING` to
compute exact near-surface distances and sweep the Eikonal equation across the
rest, which is roughly twice as fast at high resolution. Set `narrow_band = true`
and `sparse = true` to keep only samples near the surface in 4³ bricks; use
`geomlib_sdf_is_defined` to test a sample and `geomlib_sdf_make_dense` to
materialize the whole grid. Queries outside a sparse band return `FLT_MAX`.
Sparse grids save as SDF2 while dense grids keep the SDF1 layout, and
`geomlib_sdf_load` reads both.

To preserve creases, call
`geomlib_sdf_extract_surface_features(&sdf, &reference, 45.0f, true, &mesh)`.
Vertices within 1.5 voxels of a detected reference sharp edge snap onto it, and
remaining vertices optionally project onto the reference surface. Plain
`geomlib_sdf_extract_surface` rounds creases.

```c
geomlib_sdf_params params = geomlib_sdf_params_default();
params.resolution = 48;
params.algorithm = GEOMLIB_SDF_ALGORITHM_FAST_SWEEPING;
params.narrow_band = true;
params.sparse = true;
geomlib_sdf band = geomlib_mesh_to_sdf(&mesh, &params);
/* band.values is NULL; geomlib_sdf_is_defined marks the stored samples. */
geomlib_trimesh rebuilt = {0};
geomlib_sdf_extract_surface_features(&band, &mesh, 45.0f, true, &rebuilt);
geomlib_trimesh_free(&rebuilt);
geomlib_sdf_free(&band);
```

## Geometry processing

Geodesic distances use the heat method (`geomlib_compute_geodesic_from_vertex`):
a cotangent Laplacian heat solve, a normalized gradient, and a Poisson solve.
Passing `time_step <= 0` selects `mean_edge_length^2`. The public API also
provides all-pairs distances.

Mesh smoothing (`geomlib_smooth_mesh`) offers Laplacian, a two-pass Taubin
lambda/mu pair for volume preservation, and HC with the full feature-correction
term. Set `preserve_boundaries` to keep the rim fixed and `preserve_features`
to freeze vertices incident to sharp edges (`feature_angle`).

Subdivision supports Loop and modified Butterfly. Butterfly is interpolating:
regular edges use the 8-point stencil, extraordinary vertices use Zorin's masks,
and boundary edges are linear. Loop uses the cubic B-spline boundary rule
(3/4 v + 1/8(n1+n2)) and midpoint boundary edges; `preserve_boundaries` freezes
boundary vertices, and `project_to_limit` evaluates the Loop limit surface after
subdivision. `geomlib_loop_limit_position` and `geomlib_loop_limit_normal`
evaluate the limit surface with the published evaluation and tangent masks.

Surface reconstruction (`geomlib_reconstruct_surface`) selects:
- `GEOMLIB_SURFACE_GREEDY_PROJECTION`: tangent-plane one-ring fan with edge and
  normal-angle checks.
- `GEOMLIB_SURFACE_BALL_PIVOTING`: accepts triangles whose radius-`ball_radius`
  ball is empty.
- `GEOMLIB_SURFACE_MARCHING_CUBES`: oriented-point implicit field extracted with
  the SDF marching-cubes path.
- `GEOMLIB_SURFACE_ALPHA_SHAPE`: Delaunay tetrahedra with circumradius below
  `alpha`, emitting their boundary.

Ball pivoting and greedy projection use supplied normals when present and
estimate local normals otherwise.

`geomlib_decimate_qem` performs quadric-error edge collapses with link-condition
and manifold checks, optional boundary preservation, and a maximum-error cutoff.
Overdraw analysis rasterizes the mesh orthographically toward its center, culls
back faces, and reports average and maximum per-pixel fragment counts;
`geomlib_optimize_overdraw` reorders triangles front-to-back. Auto reconstruction
honors `target_faces` by choosing a marching-cubes voxel size from the point
spacing and requested face count.
Hole classification fits a PCA plane through the rim and reports planarity in
`[0,1]`. Principal curvature directions come from a height-based
normal-curvature tensor accumulated over unique incident edges. Atlas statistics
report average and maximum normalized UV edge stretch.

Collision helpers include an incremental 3D convex hull (initial tetrahedron plus
horizon expansion), OBB/AABB proxies, and QEM-decimated simplified proxies.

## CFD and tetrahedral meshes

Surface CFD functions use scalar or vector values at mesh vertices. The tetrahedral
gradient takes one scalar value per vertex and returns one gradient per cell.
Streamline tracing interpolates the velocity field barycentrically inside the
containing triangle, falling back to the nearest vertex off-surface.
Use `geomlib_fvm_build` for outward boundary faces, shared interior faces, and
cell volumes and centroids. Free its result with `geomlib_fvm_free` while keeping
the source tetrahedral mesh available for any remaining operations.

Surface tetrahedralization requires a closed, consistently oriented indexed
surface. It uses a fan when a common interior point sees the entire shell.
Otherwise it accepts a Delaunay candidate that reproduces the boundary surface
and enclosed volume; when a constrained facet is missing, Steiner points are
inserted so the boundary may be a subdivision of the input triangles while the
surface area and volume stay fixed. The constrained entry point also accepts
interior points. Use `geomlib_tetramesh_refine` to bisect eligible edges,
`geomlib_tetramesh_remove_slivers` for 2-to-3 and 3-to-2 interior flips, and
`geomlib_tetramesh_smooth` to relax interior vertices with the boundary held
fixed. The exports produce ASCII VTK legacy unstructured grids and Gmsh 2.2 meshes.

## Foundation modules (Wave 0)

New C11 modules underpin the advanced algorithms:

- `geomlib_predicates.h` — exact robust `orient2d`, `orient3d`, `incircle`,
  `insphere` and generic small-determinant signs. A floating-point filter
  falls back to expansion arithmetic; no `long double` or compiler extensions,
  so the code is identical in native and WebAssembly builds.
- `geomlib_sparse.h` — CSR matrices, Jacobi-preconditioned CG, CGLS least
  squares, dense Householder QR, symmetric Jacobi eigensolver, and the
  cotangent/mass matrices used by parameterization and spectral methods.
- `geomlib_io.h` — memory-buffer mesh I/O for OBJ, OFF, ASCII PLY and STL.
  There is no `FILE*` dependency, so the same API works in the browser.
- `geomlib_distance.h` — bounded-error symmetric Hausdorff distance with a
  rigorous upper bound derived from the 1-Lipschitz distance property.

Remeshing now projects vertices onto the original surface via a BVH closest
point query instead of scanning the first 100 original vertices.

## Wave 1: parameterization, atlas and repair

- `geomlib_parameterize_harmonic` now solves the cotangent/Tutte system with
  the sparse CG solver; `geomlib_parameterize_lscm_free` solves a real
  least-squares formulation of conformal maps with only two pinned vertices;
  `geomlib_parameterize_arap` implements the local/global ARAP solver.
  Distortion statistics now include angle RMS and L2/L∞ stretch.
- Atlas generation fills every advertised mode: planarity/distortion chart
  segmentation, MaxRects and guillotine packing, chart rotation optimization,
  and a real coverage percentage.
- Catmull–Clark subdivision is implemented and available through
  `geomlib_subdivide_mesh(..., GEOMLIB_SUBDIV_CATMULL_CLARK, ...)`.
- `geomlib_soup_to_mesh` and `geomlib_repair_nonmanifold_edges` handle
  unindexed/duplicated input and over-used edges.
- All three vertex-cache algorithms (Forsyth, Tipsify, simple greedy) work.

## Wave 2: booleans and registration

- `geomlib_mesh_boolean` performs union, intersection and difference with a
  shared-grid SDF min/max combination and marching-cubes extraction.
  `voxel_hint <= 0` chooses one grid cell per 1/32 of the combined bounds.
- `geomlib_icp_rigid` performs rigid point-to-surface ICP and returns the final
  RMS and a column-major 4x4 transform; `geomlib_transform_mesh` applies it.
- `geomlib_blue_noise_surface` generates Poisson-disk samples on a surface with
  area-weighted triangle selection and deterministic RNG support.
- `geomlib_reconstruct_poisson` reconstructs a surface from oriented points by
  solving the discrete Poisson equation for an indicator field.
- `geomlib_remesh_acvd` performs CVT/Lloyd relaxation of surface sites and
  retriangulates them, producing a near-uniform mesh with a target vertex count.
- `geomlib_deform_arap` deforms a surface with position handles using the local/
  global as-rigid-as-possible solver.
- `geomlib_compute_manifold_harmonics` returns the lowest Laplace-Beltrami
  eigenvalues and eigenfunctions of a mesh.
- `geomlib_detect_shapes_ransac` finds planes, spheres and cylinders in point
  clouds with deterministic seeded RANSAC.
- `geomlib_remesh_to_quads` converts triangle meshes into mixed quad/triangle
  meshes with dihedral and angle quality constraints.
- `geomlib_remesh_anisotropic` refines meshes adaptively using curvature-derived
  metric tensors.
- `geomlib_sdf_marching_tets` extracts watertight isosurfaces via six tetrahedra
  per grid cube (no ambiguous faces).
- `geomlib_sdf_dual_contouring` extracts isosurfaces with QEM-based vertex
  placement, preserving sharp features.
- `geomlib_compute_saliency` and `geomlib_segment_by_saliency` measure and
  partition shape "interestingness" from multi-scale curvature.
- `geomlib_extract_curve_skeleton` contracts a surface into a 1D topology-
  preserving curve skeleton; `geomlib_compute_straight_skeleton` computes the
  2D straight skeleton of a simple polygon; `geomlib_offset_polygon` offsets
  simple 2D polygons inward or outward; `geomlib_vsa_simplify` approximates a
  mesh with planar proxies and returns a coarse simplified mesh.
- `geomlib_pointcloud_bilateral_smooth` and
  `geomlib_pointcloud_remove_outliers` denoise and clean point clouds while
  preserving features; `geomlib_alpha_wrap` wraps arbitrary geometry in a
  watertight offset surface; `geomlib_slice_mesh` intersects a mesh with a
  plane and returns polylines; `geomlib_approx_convex_decomposition` splits a
  mesh into convex hull pieces; `geomlib_fair_mesh` performs unconditionally
  stable implicit mean-curvature flow smoothing; `geomlib_fast_winding_number`
  accelerates generalized winding-number queries with a cluster hierarchy.
  Binary little-endian PLY is supported via the `"plyb"` format name;
  `geomlib_partition_mesh` creates balanced, spatially coherent vertex parts;
  `geomlib_decimate_vertex_clustering` provides a fast lossy simplification;
  `geomlib_marching_squares` extracts 2D zero contours from scalar grids;
  `geomlib_delaunay2d` triangulates 2D point sets;
  `geomlib_geodesic_fast_marching` computes surface distance fields;
  `geomlib_geodesic_length_unfolding` gives continuous path lengths;
  `geomlib_skin_mesh` and `geomlib_skin_mesh_dualquat` apply linear-blend and
  dual-quaternion skinning transforms; `geomlib_compute_cross_field` generates
  combed 4-directional fields for quad/hex meshing;
  `geomlib_compute_shape_diameter` measures local thickness for segmentation;
  `geomlib_kdtree_*` accelerates point-set kNN and radius queries;
  `geomlib_bvh_refit` updates a BVH after vertex deformation;
  `geomlib_octree_*` accelerates point-set radius and kNN queries;
  `geomlib_solve_sparse_direct` solves sparse SPD systems directly;
  `geomlib_hex_dominant` merges tetrahedra into hex-dominant volume meshes;
  `geomlib_reconstruct_rbf` reconstructs surfaces with compactly supported RBFs;
  `geomlib_fluid_*` provides grid-based incompressible fluid advection/projection;
  `geomlib_compute_voronoi2d` builds clipped 2D Voronoi cells with nearest-site queries;
  `geomlib_voxelize_mesh` rasterizes triangle meshes into surface and volume voxels;
  `geomlib_spatial_sort2/3` provide Morton and Hilbert cache-coherent point orderings;
  `geomlib_mesh_mesh_intersection` extracts cross-section segments between two meshes;
  `geomlib_cdt2d_build` triangulates 2D domains with preserved constraint segments;
  `geomlib_offset_mesh` grows or shrinks closed meshes via signed distance levels;
  `geomlib_farthest_point_sampling` spreads samples over a mesh with coverage radii;
  `geomlib_medial_axis2d` extracts a polygon's medial axis from boundary-sample Voronoi cells;
  `geomlib_polygon_boolean` computes 2D polygon union, intersection and difference;
  `geomlib_nonrigid_icp` fits smooth deformations between triangle meshes;
  `geomlib_convex_hull2d` computes planar convex hulls with monotone chain;
  `geomlib_minkowski_sum2d` computes convex Minkowski sums of planar shapes;
  `geomlib_triangulate_with_holes` triangulates polygons containing holes;
  `geomlib_simplify_polyline` reduces polylines with Douglas-Peucker;
  `geomlib_mls_project` projects queries onto point-cloud MLS surfaces;
  `geomlib_sdf_raycast` sphere-traces rays through signed distance fields;
  `geomlib_geodesic_dijkstra` gives exact mesh-edge shortest-path distances;
  `geomlib_parameterize_abf` reduces 3D-to-2D corner-angle distortion;
  `geomlib_deform_mean_value2d` warps points with affine-precise cage coordinates;
  `geomlib_mesh_bilateral_smooth` denoises meshes with normal-weighted averaging;
  `geomlib_resample_catmull_rom` resamples curves with uniform or centripetal splines;
  `geomlib_improve_mesh_quality` flips edges to improve triangle angles;
  `geomlib_refine_edges` subdivides long edges for isotropic refinement;
  `geomlib_coarsen_edges` collapses short edges with manifold checks;
  `geomlib_refine_polygon2d` refines planar domains with Delaunay Steiner points;
  `geomlib_surface_voronoi_compute` partitions meshes into geodesic Voronoi regions;
  `geomlib_icp_point_to_plane` aligns meshes with point-to-plane ICP;
  `geomlib_mesh_contour_lines` extracts scalar-field isolines on meshes;
  `geomlib_mesh_region_growing` segments faces by dihedral continuity;
  `geomlib_pointcloud_curvature` estimates local curvature and normals by PCA;
  `geomlib_compute_hks` computes multi-scale heat kernel shape descriptors;
  `geomlib_convex_decompose2d` splits simple polygons into convex parts;
  `geomlib_bounding_sphere_points/mesh` compute Ritter bounding spheres;
  `geomlib_minimum_enclosing_sphere` computes the exact smallest enclosing ball;
  `geomlib_compute_mass_properties` integrates volume, centroid and inertia;
  `geomlib_largest_inscribed_circle` finds a polygon's Chebyshev centre;
  `geomlib_minimum_area_rectangle` computes the smallest enclosing rectangle;
  `geomlib_polyline_mesh_intersection` finds surface crossings of 3D polylines;
  `geomlib_mesh_ambient_occlusion` estimates vertex visibility by hemisphere rays;
  `geomlib_spectral_segment` partitions meshes from their spectral embedding;
  `geomlib_mesh_mesh_distance` finds closest points between two surfaces;
  `geomlib_detect_symmetry` finds principal reflection planes of a mesh;
  `geomlib_pointcloud_kmeans` clusters point sets with deterministic Lloyd iterations;
  `geomlib_pointcloud_dbscan` finds density clusters and marks noise;
  `geomlib_resample_polyline` redistributes curve points by arc length;
  `geomlib_compute_polygon_moments` integrates planar area, centroid and inertia;
  `geomlib_clip_polygon_convex` clips polygons against convex regions;
  `geomlib_minimum_enclosing_circle` computes exact smallest enclosing circles;
  `geomlib_geodesic_path_between` extracts shortest-path vertex chains on meshes;
  `geomlib_extract_geodesic_ball` selects vertices within a geodesic radius;
  `geomlib_fit_rigid_points` fits rigid/similarity transforms between point sets;
  `geomlib_fit_circle` fits circles to planar point sets;
  `geomlib_fit_sphere` fits spheres to 3D point sets;
  `geomlib_fit_plane` fits planes by PCA;
  `geomlib_chaikin_subdivide` smooths polylines by corner cutting;
  `geomlib_fit_ellipse` fits ellipses to planar point sets;
  `geomlib_fit_cylinder` fits infinite cylinders to 3D points;
  `geomlib_segment_closest_points` computes closest points between 3D segments;
  `geomlib_ray_plane_intersection` intersects rays and segments with planes;
  `geomlib_triangle_intersection` finds intersection segments of triangles;
  `geomlib_polygon_boolean_holes` applies Boolean ops to polygons with holes;
  `geomlib_geodesic_diameter_exact/sweep` find mesh diameter vertices;
  `geomlib_deform_tetra_cage` warps points with tetrahedral barycentric weights;
  `geomlib_poisson_disk_volume` generates blue-noise 3D point sets;
  `geomlib_point_in_polygon` provides containment, winding and signed distance;
  `geomlib_build_laplacian` assembles cotangent stiffness and lumped mass;
  `geomlib_quat_*` provides quaternion algebra, slerp and rotation;
  `geomlib_clip_mesh_plane` cuts meshes by planes and caps the section;
  `geomlib_segment2d_intersect` classifies planar segment intersections;
  `geomlib_compute_mst` builds minimum spanning trees of point sets;
  `geomlib_compute_pca` finds principal axes of point clouds;
  `geomlib_sphere_sphere_intersection` computes intersection circles;
  `geomlib_frechet_distance` measures discrete curve similarity;
  `geomlib_poisson_disk_polygon` samples blue-noise points inside polygons;
  `geomlib_smooth_polyline` smooths curves by Laplacian flow;
  `geomlib_cage2d_harmonic_*` deforms points with harmonic cage coordinates;
  `geomlib_rasterize_polygon2d` samples polygons onto grids with signed distance;
  `geomlib_visibility_shortest_path` plans paths around polygon obstacles;
  `geomlib_aabb_overlap_pairs` finds broad-phase overlaps between meshes;
  `geomlib_solid_angle_mesh` computes exact solid angle and winding numbers;
  `geomlib_simplify_polyline_vw` simplifies curves by effective area;
  `geomlib_bezier_evaluate` evaluates and subdivides Bezier curves;
  `geomlib_ray_sphere_intersect` computes analytic ray-sphere roots;
  `geomlib_plane_plane_intersection` finds the line where two planes meet;
  `geomlib_sphere_aabb_overlap` tests sphere-box overlap with penetration depth;
  `geomlib_circle2d_intersect` computes planar circle intersections;
  `geomlib_segment_triangle_intersect` intersects finite segments with triangles;
  `geomlib_ray_cylinder_intersect` intersects rays with infinite cylinders;
  `geomlib_capsule_capsule_overlap` tests capsule collision with penetration depth;
  `geomlib_sphere_triangle_overlap` tests sphere-triangle penetration;
  `geomlib_triangle_aabb_overlap` tests triangle-box overlap via SAT;
  `geomlib_convex_polygon_diameter` finds farthest convex-polygon vertices;
  `geomlib_convex_polygon_min_width` computes the minimum supporting-line width;
  `geomlib_aabb_distance` finds the gap and closest points between boxes;
  `geomlib_ray_obb_intersect` intersects rays with oriented bounding boxes;
  `geomlib_triangle_triangle_distance` finds closest points between triangles;
  `geomlib_plane_sphere_intersect` computes plane-sphere intersection circles;
  `geomlib_circle_line_intersect` intersects planar circles with lines and segments;
  `geomlib_ray_tetrahedron_intersect` intersects rays with tetrahedra;
  `geomlib_polyline_closest_point` projects points onto polylines;
  `geomlib_capsule_aabb_overlap` tests capsules against boxes;
  `geomlib_polygon_distance2d` finds closest features of planar polygons;
  `geomlib_triangle_tetrahedron_intersect` tests triangle-solid overlap;
  `geomlib_sphere_capsule_overlap` tests spheres against capsules;
  `geomlib_tetrahedra_intersect` tests tetrahedron-tetrahedron overlap;
  `geomlib_sphere_obb_overlap` tests spheres against oriented boxes;
  `geomlib_circle_aabb_overlap` tests circles against planar rectangles;
  `geomlib_capsule_obb_overlap` tests capsules against oriented boxes;
  `geomlib_ray_cone_intersect` intersects rays with cones;
  `geomlib_triangle_obb_overlap` tests triangles against oriented boxes;
  `geomlib_aabb_obb_overlap` tests axis-aligned against oriented boxes;
  `geomlib_obb_obb_overlap` tests two oriented boxes with SAT;
  `geomlib_point_in_ellipse` queries oriented 2D ellipses;
  `geomlib_capsule_triangle_overlap` tests capsules against triangles;
  `geomlib_ellipse_closest_point` projects points onto ellipse boundaries;
  `geomlib_ellipse_perimeter` approximates ellipse perimeter via Ramanujan;
  `geomlib_circle_triangle_intersect` tests planar circles against triangles;
  `geomlib_line2d_intersect` computes planar line intersections;
  `geomlib_polyline_arc_table` builds cumulative arc-length parameterizations;
  `geomlib_circle_intersection_area` computes circle lens and union areas;
  `geomlib_polygon_intersection_area` computes convex polygon overlap area;
  `geomlib_circle_ellipse_overlap` tests circles against oriented ellipses;

## Clean-room process and Lean formalization

Algorithms are implemented from papers only. Every new module has a contract
in [`docs/contracts/`](docs/contracts/README.md) recording its source,
specification, derivation, API, tests and formalization target.

The [`formal/`](formal/) directory contains core-Lean (mathlib-free) machine
checked specifications:

- `GeomLibFormal/MarchingCubes.lean` — verified cube-face parity, crossing
  parity, complement symmetry and bound; `formal/gen_mc_table.lean` emits
  `include/geomlib/generated/mc_table.h` from those definitions.
- `GeomLibFormal/Predicates.lean` — exact integer determinant specifications.
- `GeomLibFormal/Delaunay.lean` — in-circle translation, swap, scaling and
  rotation axioms verified over a bounded integer domain.
- `GeomLibFormal/Topology.lean` — combinatorial map inverse and ∂² cancellation.

Check them with `lean formal/GeomLibFormal/*.lean` (Lean 4.36). A
`lakefile.lean` declaring mathlib is provided for larger future developments,
but mathlib is not required by the checked files.

## Building WebAssembly

```sh
./scripts/configure.sh --wasm
cd web/build && emmake make
cd ../.. && node web/smoke.mjs
```

Both native and WASM builds share the same source list
(`cmake/geomlib_sources.cmake`) and compile the C core as strict C11
(`CMAKE_C_EXTENSIONS OFF`). Use `--strict-c11` to add `-pedantic -Werror`.

## License

MIT License - See [LICENSE](LICENSE) for details.
