# Contract: fast-marching-style surface geodesics

**Status:** implemented (Wave 3).

## Sources (papers only)

- Sethian, "A Fast Marching Level Set Method for Monotonically Advancing
  Fronts", 1996.
- Kimmel & Sethian, "Computing Geodesic Paths on Manifolds", PNAS 1998.

## Approach

Dijkstra on the mesh edge graph gives an initial upper bound. Repeated
triangle-interior relaxation samples points on each triangle edge, interpolates
their distances from the endpoints, and relaxes the opposite vertex by the
interpolated distance plus the in-triangle Euclidean distance. This approximates
wavefront propagation across triangle interiors.

## API / constraints

`geomlib_geodesic_fast_marching(mesh, source, passes, out)`. C11, wasm32 safe.

## Acceptance

`tests/test_geodesic_fmm.cpp`: unit sphere from one vertex gives a finite
distance field with max 3.1365 (analytic diameter ≈ π). WASM smoke checks the
returned array and zero source distance.

## Lean target

Not yet formalized.
