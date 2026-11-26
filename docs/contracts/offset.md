# Contract: SDF-based mesh offsetting

**Status:** implemented (Wave 3).

## Source (paper only)

- Barentzen et al., "Fast and Robust Offsetting of Polygonal Meshes", 2015.
- Standard signed-distance level-set offsetting.

## Approach

Sample a signed distance field (negative inside) on a grid padded around the
input mesh, then extract the iso level equal to the requested offset. Positive
offsets grow the solid and negative offsets shrink it. Grid voxel size is
derived from the mesh extent and the requested resolution/padding.

## API / constraints

`geomlib_offset_mesh(mesh, offset, params, out)` with
`geomlib_offset_params_default`. Requires a closed, consistently oriented mesh.
Large negative offsets can collapse thin features. C11, wasm32 safe.

## Acceptance

`tests/test_offset.cpp`: a radius-0.5 sphere offset by +0.1 gives radius 0.6117
and volume 0.8838 (analytic 0.905); offset by -0.1 gives radius 0.4412 and
volume 0.2591 (analytic 0.268). WASM smoke checks an outward offset.

## Lean target

Not yet formalized.
