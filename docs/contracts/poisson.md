# Contract: Poisson surface reconstruction

**Status:** implemented (Wave 2).

## Source (paper only)

- Kazhdan, Bolitho, Hoppe, "Poisson Surface Reconstruction", Eurographics
  Symposium on Geometry Processing 2006.

## Approach

Voxelized Poisson reconstruction: oriented points are splatted with Gaussian
weights into a 3D vector field of normals; the divergence is computed with
clamped central differences; a 7-point Neumann Laplacian system is solved for
the indicator function with the Jacobi-preconditioned sparse CG solver (one
corner pinned to remove the constant nullspace); the iso level is the median
indicator value at the input points; marching cubes extracts the isosurface.

Unlike the existing oriented-point implicit-field reconstruction in
`geomlib_surface.c`, this actually solves the Poisson equation.

## API / constraints

`geomlib_reconstruct_poisson(cloud, params, out)`; default resolution 48,
capped at 128 per axis. Requires normals in the point cloud. C11, wasm32 safe.

## Acceptance

`tests/test_poisson.cpp`: 400 analytic unit-sphere points reconstruct to a
26000-triangle mesh with volume 4.1995 (exact sphere volume 4.1888).
WASM smoke reconstructs from 200 oriented points.

## Lean target

Discrete Poisson consistency/convergence identities are future work.
