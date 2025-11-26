# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

### Configure and Build
```bash
# Quick configure and build (uses scripts/configure.sh)
./scripts/configure.sh
cd build && make -j8

# Manual configuration
mkdir build && cd build
cmake ..
make -j8
```

### Build Options
- `GEOMLIB_BUILD_TESTS=ON|OFF` - Build tests (default: ON)
- `GEOMLIB_BUILD_EXAMPLES=ON|OFF` - Build examples (default: OFF)

### Configure Script Usage
```bash
# Debug build
./scripts/configure.sh --debug

# Custom build directory
./scripts/configure.sh --build-dir build-debug

# With examples
./scripts/configure.sh --examples

# No tests
./scripts/configure.sh --no-tests

# Custom install prefix
./scripts/configure.sh --prefix /usr/local
```

## Testing

### Run All Tests
```bash
cd build && ctest
# or
cd build && ctest --verbose
```

### Run Individual Tests
```bash
cd build/tests
./test_linalg      # Linear algebra tests
./test_solver      # Linear solver tests
./test_mesh        # Mesh processing tests
./test_random      # Random number generation tests
./test_montecarlo  # Monte Carlo geometry processing tests
./test_bvh         # BVH and ray tracing tests
```

## Architecture

### Core Design Principles
- **C++17 compatible** but core implementations in C for maximum portability
- **No RTTI, no C++ exceptions** - suitable for embedded and real-time systems
- **No external dependencies** - custom linear algebra (no Eigen or similar libraries)
- **Column-major matrices** - following OpenGL convention

### Six-Module Architecture

**1. Linear Algebra (`geomlib_linalg.h` / `geomlib_linalg.c`)**
- Fundamental data types: `geomlib_vec2/3/4`, `geomlib_mat3/4`
- Double precision variants: `geomlib_vec2d/3d/4d`
- Matrix operations use column-major order
- Transformation matrices: translate, scale, rotate, look_at, perspective, ortho
- All operations are pure C functions with `extern "C"` linkage

**2. Linear Solvers (`geomlib_solver.h` / `geomlib_solver.c`)**
- Dense solvers: Gaussian elimination, LU decomposition
- Iterative solvers: Conjugate Gradient (for symmetric positive-definite), Jacobi
- Eigenvalue: Power iteration for dominant eigenvalue
- Least squares solver
- Return codes: `GEOMLIB_SOLVER_SUCCESS`, `GEOMLIB_SOLVER_ERROR_SINGULAR`, etc.

**3. Mesh Processing (`geomlib_mesh.h` / `geomlib_mesh.c`)**
- Data structures: `geomlib_polygon2d/3d`, `geomlib_trimesh`, `geomlib_polymesh`, `geomlib_aabb`
- Polygon operations: area, perimeter, centroid, convexity testing, point containment
- Mesh operations: add vertices/triangles, compute normals, compute bounds
- Primitive generation: `geomlib_create_box/sphere/cylinder/plane`
- Mesh conversion: triangulate polygon meshes

**4. Random Number Generation (`geomlib_random.h` / `geomlib_random.c`)**
- PCG32 (Permuted Congruential Generator) - fast, statistically good PRNG
- State: `geomlib_pcg32` with seed and sequence initialization
- Integer generation: `geomlib_pcg32_next()` for uint32_t
- Float/double generation: `geomlib_pcg32_nextf/d()` using IEEE 754 bit manipulation (no division)
- Range functions: `geomlib_pcg32_range_int/float/double()` for bounded values
- Technique from https://github.com/wjakob/pcg32 - efficient [0,1) generation via bit patterns

**5. Monte Carlo Geometry Processing (`geomlib_montecarlo.h` / `geomlib_montecarlo.c`)**
- **Grid-free approach** to geometry processing avoiding mesh discretization
- **Triangle sampling**: Uniform sampling using barycentric coordinates (standard method from pbr-book.org)
- **Mesh surface sampling**: Area-weighted sampling with CDF for O(log n) selection
- **Monte Carlo integration**: Estimators for surface integrals with variance computation
- **Closest point queries**: Triangle and mesh closest point (foundation for Walk on Spheres)
- **Utility samplers**: Sphere surface/volume, hemisphere (uniform and cosine-weighted)
- Based on [Keenan Crane's Monte Carlo Geometry Processing](https://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/)

**6. BVH and Ray Tracing (`geomlib_bvh.h` / `geomlib_bvh.c`)**
- **Multi-level BVH** with Surface Area Heuristic (SAH) construction
- **Ray-triangle intersection**: Möller-Trumbore algorithm (watertight, no plane equation)
- **Ray-AABB intersection**: Optimized slab method for traversal
- **Ray queries**: Closest hit, occlusion (any-hit), normal interpolation
- **Point queries**: Closest point, k-NN, radius search using BVH acceleration
- **Quality levels**: Fast/Low/Medium/High build modes (SAH buckets: 12, max leaf primitives: 4)
- Based on [pbr-book.org BVH chapter](https://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies)

### Memory Management Pattern
The library uses manual memory management:
- Mesh structures have explicit `_init()` and `_free()` functions
- Example: `geomlib_trimesh_init()` / `geomlib_trimesh_free()`
- Always pair init/free calls to prevent memory leaks

### Test Structure
Tests are simple C++ executables (no external test framework):
- `test_linalg.cpp` - Vector, matrix, transformation tests
- `test_solver.cpp` - Linear solver correctness tests
- `test_mesh.cpp` - Polygon, mesh, primitive generation tests
- `test_random.cpp` - PCG32 PRNG tests (determinism, distribution, ranges)
- `test_montecarlo.cpp` - Monte Carlo sampling and integration tests
- `test_bvh.cpp` - BVH construction, ray tracing, point query tests
- Tests use basic assertions and return 0 for success

## Code Organization

```
geomlib/
├── include/geomlib/       # Public headers
│   ├── geomlib.h          # Main header (includes all modules)
│   ├── geomlib_linalg.h      # Linear algebra types and operations
│   ├── geomlib_solver.h      # Linear system solvers
│   ├── geomlib_mesh.h        # Mesh/polygon data structures
│   ├── geomlib_random.h      # PCG32 random number generation
│   ├── geomlib_montecarlo.h  # Monte Carlo geometry processing
│   └── geomlib_bvh.h         # BVH and ray tracing
├── src/                      # C implementations
│   ├── geomlib_linalg.c      # ~600 lines
│   ├── geomlib_solver.c      # ~500 lines
│   ├── geomlib_mesh.c        # ~900 lines
│   ├── geomlib_random.c      # ~150 lines
│   ├── geomlib_montecarlo.c  # ~600 lines
│   └── geomlib_bvh.c         # ~800 lines
├── tests/                 # C++ test files
├── scripts/               # Build scripts
│   └── configure.sh       # CMake configuration helper
└── build/                 # Build directory (gitignored)
```

## Adding New Features

### Adding Vector/Matrix Operations
- Add declaration to `include/geomlib/geomlib_linalg.h`
- Implement in C in `src/geomlib_linalg.c`
- Add test cases to `tests/test_linalg.cpp`
- Use `extern "C"` for C++ compatibility

### Adding New Solvers
- Add declaration to `include/geomlib/geomlib_solver.h` with appropriate return codes
- Implement in C in `src/geomlib_solver.c`
- Add test cases to `tests/test_solver.cpp`
- Document convergence criteria and limitations

### Adding Mesh Operations
- Add declaration to `include/geomlib/geomlib_mesh.h`
- Implement in C in `src/geomlib_mesh.c`
- Follow existing memory management patterns (init/free)
- Add test cases to `tests/test_mesh.cpp`

### Adding Random Number Features
- Add declaration to `include/geomlib/geomlib_random.h`
- Implement in C in `src/geomlib_random.c`
- PCG32 state is lightweight (`geomlib_pcg32` struct with state and inc)
- Add test cases to `tests/test_random.cpp` with statistical validation

### Adding Monte Carlo Features
- Add declaration to `include/geomlib/geomlib_montecarlo.h`
- Implement in C in `src/geomlib_montecarlo.c`
- Sampling functions require PCG32 RNG state
- Mesh sampling requires preprocessing with `geomlib_mesh_build_sample_distribution()`
- Add test cases to `tests/test_montecarlo.cpp` with statistical validation
- Key patterns: uniform sampling uses barycentric coordinates, area-weighted uses CDF

## Monte Carlo Geometry Processing

### Core Concepts
This module implements foundational Monte Carlo techniques for geometry processing, inspired by the grid-free approach that avoids spatial discretization.

### Triangle Sampling
- **Uniform sampling**: `geomlib_sample_triangle_uniform()` uses the standard barycentric method
- **Formula**: `b0 = 1 - sqrt(u), b1 = v * sqrt(u), b2 = 1 - b0 - b1`
- **Why it works**: The sqrt transformation prevents clustering near vertices
- **PDF**: For uniform sampling, `pdf = 1 / triangle_area`

### Mesh Surface Sampling
1. **Preprocessing**: Build CDF with `geomlib_mesh_build_sample_distribution()`
   - Computes area for each triangle
   - Builds cumulative distribution function (CDF)
   - O(n) preprocessing, O(log n) sampling via binary search
2. **Sampling**: `geomlib_sample_mesh_surface()` samples uniformly across entire surface
3. **Memory**: Remember to free distribution with `geomlib_mesh_free_sample_distribution()`

### Monte Carlo Integration
- **Estimator**: `∫_S f(x) dA ≈ (1/N) * Σ f(x_i) / pdf(x_i)`
- **Uniform case**: For uniform sampling, simplifies to `(area/N) * Σ f(x_i)`
- **Variance**: Use `geomlib_mc_estimate_variance()` to estimate error
- **Example**: Estimate π by sampling unit square and counting points in unit circle

### Closest Point Queries
- **Triangle**: `geomlib_closest_point_triangle()` handles all regions (interior, edges, vertices)
- **Mesh**: `geomlib_closest_point_mesh()` brute-force search (consider BVH for production)
- **Purpose**: Foundation for Walk on Spheres algorithm (grid-free PDE solving)

### Utility Samplers
- **Sphere surface**: `geomlib_sample_sphere_uniform()` - spherical coordinates method
- **Sphere volume**: `geomlib_sample_sphere_volume()` - rejection sampling
- **Hemisphere**: `geomlib_sample_hemisphere_uniform()` - for geometric queries
- **Cosine hemisphere**: `geomlib_sample_hemisphere_cosine()` - for diffuse reflection (rendering)

## BVH and Ray Tracing

### Core Concepts
The BVH (Bounding Volume Hierarchy) module provides efficient spatial acceleration for ray tracing and point queries using a binary tree of axis-aligned bounding boxes.

### BVH Construction
1. **Build Quality Settings**:
   - `GEOMLIB_BVH_QUALITY_FAST` - Median split, fastest build
   - `GEOMLIB_BVH_QUALITY_LOW` - Basic SAH with fewer buckets
   - `GEOMLIB_BVH_QUALITY_MEDIUM` - Balanced SAH (12 buckets, recommended)
   - `GEOMLIB_BVH_QUALITY_HIGH` - Full SAH with more analysis

2. **SAH (Surface Area Heuristic)**:
   - Evaluates split quality using surface area as proxy for ray intersection probability
   - Cost model: `traversal_cost + (left_SA/total_SA * left_count + right_SA/total_SA * right_count) * intersection_cost`
   - Uses 12 buckets along chosen axis for O(n) split evaluation
   - Falls back to median split for small primitive counts (≤4)

3. **Build Process**:
   ```c
   geomlib_bvh bvh;
   geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
   // Use BVH...
   geomlib_bvh_free(&bvh);
   ```

### Ray Tracing
- **Möller-Trumbore intersection**: Fast watertight ray-triangle test without plane equation
- **Stack-based traversal**: 64-level stack for iterative BVH traversal (no recursion)
- **Closest hit**: `geomlib_bvh_trace_ray()` finds nearest intersection
- **Occlusion**: `geomlib_bvh_test_occlusion()` early-exits on first hit (shadow rays)
- **Normal interpolation**: Automatically interpolates vertex normals if available

### Point Queries
- **Closest point**: `geomlib_bvh_closest_point()` finds nearest surface point
  - Uses BVH to prune search space
  - Much faster than brute-force for large meshes
- **k-NN**: `geomlib_bvh_knn_query()` finds k nearest neighbors (basic implementation)
- **Radius search**: `geomlib_bvh_radius_query()` finds all points within radius
  - Culls entire BVH nodes outside search sphere
  - Returns up to max_results points

### Performance Characteristics
- **Build time**: O(n log n) with SAH, O(n log n) with median
- **Memory**: ~2n nodes for n triangles (binary tree)
- **Ray traversal**: O(log n) average case
- **Point queries**: O(log n) with BVH vs O(n) brute-force

### Usage Example
```c
// Build BVH
geomlib_bvh bvh;
geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);

// Ray tracing
geomlib_ray ray = geomlib_ray_create(origin, direction);
geomlib_ray_hit hit;
if (geomlib_bvh_trace_ray(&bvh, &ray, &hit)) {
    printf("Hit at t=%.2f, normal=(%.2f, %.2f, %.2f)\n",
           hit.t, hit.normal.x, hit.normal.y, hit.normal.z);
}

// Closest point query
geomlib_vec3 closest;
float distance;
geomlib_bvh_closest_point(&bvh, query, &closest, &distance, NULL);

geomlib_bvh_free(&bvh);
```

## Compiler Flags

The build system automatically sets:
- `-fno-rtti -fno-exceptions` (GCC/Clang)
- `/GR- /EHsc` (MSVC)
- C11 standard for C files
- C++17 standard for C++ files
- Links with `-lm` on Unix platforms
