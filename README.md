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

## License

MIT License - See [LICENSE](LICENSE) for details.
