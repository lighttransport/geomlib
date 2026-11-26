# Contract: variational RBF surface reconstruction

**Status:** implemented (Wave 3).

## Sources (papers only)

- Carr, Snoeyink, Axen, "Computing and Rendering Point Set Surfaces", 2001.
- Turk & O'Brien, "Shape Transformation Using Variational Implicit Functions",
  1999.

## Approach

Fit compactly supported Wendland C2 radial basis functions to input points by
solving a dense linear system with unit right-hand side, evaluate the implicit
field on a dense grid, take a high iso level close to the data points, and
extract the surface with marching cubes. Normals are optional.

## API / constraints

`geomlib_reconstruct_rbf(cloud, params, out)`; resolution capped at 96 and
point count at 500 by default. C11, wasm32 safe.

## Acceptance

`tests/test_rbf_reconstruct.cpp`: 200 sphere points reconstruct to 4252
triangles with volume 4.75 (analytic 4.19). WASM smoke checks a smaller cloud.

## Lean target

Not yet formalized.
