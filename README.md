# geomlib

A portable geometry processing library in C++17 with core functions implemented in C.

## Features

- **Linear Algebra**: Vectors (vec2, vec3, vec4), matrices (mat3, mat4), transformations
- **Linear Solvers**: Gaussian elimination, LU decomposition, Conjugate Gradient, Jacobi iteration
- **Mesh Processing**: Triangle meshes, polygon meshes, primitive generation
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

## License

MIT License - See [LICENSE](LICENSE) for details.
