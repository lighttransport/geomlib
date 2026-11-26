# Contract: tetrahedral cage deformation

**Status:** implemented (Wave 3).

## Source (paper only)

- Barycentric coordinates of a tetrahedron via signed sub-volume ratios;
  affine-precise cage deformation.

## Approach

Compute barycentric weights as ratios of the four sub-tetrahedron volumes to
the cage's volume, then apply the deformed cage vertices with those weights.
Works for points inside and outside the cage (using signed volumes).

## API / constraints

`geomlib_tetra_weights(cage, q, weights)` and
`geomlib_deform_tetra_cage(cage, deformed, points, n, out)`. Requires a
non-degenerate tetrahedron. C11, wasm32 safe.

## Acceptance

`tests/test_cage3d.cpp`: weights for an interior point sum to 1, the identity
cage leaves points unchanged, and a scale-2 + translate cage reproduces the
affine image exactly. WASM smoke checks the affine case.

## Lean target

Not yet formalized.
