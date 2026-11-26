# Contract: triangle-mesh voxelization

**Status:** implemented (Wave 3).

## Source (paper only)

- Akenine-Moller, "Fast 3D Triangle-Box Overlap Testing", 2001 (13-axis SAT).
- Standard ray-parity inside test for closed surfaces.

## Approach

For each triangle, enumerate candidate voxels from its axis-aligned bounds and
keep those whose centre-box overlaps the triangle under the 13 separating axes
(3 box normals, triangle normal, 9 cross products). For closed meshes, cast a
ray along +X from every voxel centre through a BVH and count crossings; odd
means inside. The grid spans the mesh AABB with a 5% margin.

## API / constraints

`geomlib_voxelize_mesh(mesh, nx, ny, nz)` and `geomlib_voxel_grid_free`.
The surface pass is conservative; the volume pass requires a closed, consistently
oriented mesh. Total voxels capped at 64M. C11, wasm32 safe.

## Acceptance

`tests/test_voxelize.cpp`: a unit box at 24³ gives volume estimate 1.0252 and
2648 surface voxels; a radius-0.5 sphere gives 0.5030 vs analytic 0.5236.
WASM smoke checks both label arrays are produced.

## Lean target

Not yet formalized.
