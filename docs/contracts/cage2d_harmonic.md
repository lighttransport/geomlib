# Contract: 2D harmonic cage coordinates

**Status:** implemented (Wave 3).

## Source (paper only)

- Floater, "Mean Value Coordinates", CAGD 2003; harmonic-like angle/inverse-
  distance weights for convex polygons.

## Approach

Compute signed angles subtended at the query point by each incident pair of
cage edges, take tangent half-angles, divide by the distance to the cage vertex,
and normalize. Deformation applies the affine combination of moved cage
vertices.

## API / constraints

`geomlib_cage2d_harmonic_weights(cage, n, q, weights)` and
`geomlib_cage2d_harmonic_deform(...)`. Affine precision holds for convex cages.
C11, wasm32 safe.

## Acceptance

`tests/test_cage2d_harmonic.cpp`: weights sum to one, identity cage is a no-op,
and a scale/translate cage deformation is reproduced exactly. WASM smoke checks
a translation.

## Lean target

Not yet formalized.
