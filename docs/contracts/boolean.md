# Contract: mesh booleans

**Status:** implemented (Wave 2).

## Approach

Grid-based signed-distance booleans. Both solids are sampled onto one shared
dense grid; union is `min(dA,dB)`, intersection is `max(dA,dB)`, difference
`A-B` is `max(dA,-dB)`; the result is extracted with the existing marching
cubes path. This is approximate (voxel_hint controls accuracy) but robust for
touching, nested and self-intersecting inputs.

A first surface-fragment implementation was attempted and rejected: winding
numbers at coincident fragments are ~0.5 and fragment stitching produced
non-watertight difference results. The grid method is the retained approach.

## API / constraints

`geomlib_mesh_boolean(a, b, op, voxel_hint, out)`. C11, uses BVH and SDF
modules; `voxel_hint <= 0` selects max_extent / 32. WASM: yes (resolution
capped at 256 per axis).

## Acceptance

`tests/test_boolean.cpp`: two unit boxes offset by 0.5 produce volumes
1.5 / 0.5 / 0.5 for union / intersection / difference within 0.01.
WASM smoke checks a box-sphere union.

## Lean target

Set semantics of booleans (`solid(A ∩ B) = solid(A) ∩ solid(B)`) is a
research-grade future formalization.
