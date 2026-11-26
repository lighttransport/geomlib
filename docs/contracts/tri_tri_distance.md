# Contract: triangle–triangle distance

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-feature enumeration between two triangles; see Ericson,
  "Real-Time Collision Detection".

## Approach

Return 0 when the triangles intersect; otherwise enumerate all vertex-to-
triangle projections (6) and all edge-to-edge segment closest pairs (9), then
keep the smallest distance with its witnesses.

## API / constraints

`geomlib_triangle_triangle_distance(a0,a1,a2,b0,b1,b2,out)` fills distance and
closest points. C11, wasm32 safe.

## Acceptance

`tests/test_tri_tri_distance.cpp`: triangles 3 apart along x give distance 2 at
vertices (1,0,0) and (3,0,0); overlapping triangles give 0; triangles in
parallel planes two apart give 2. WASM smoke checks the separated case.

## Note

The shared triangle intersection test's separation check was corrected to test
each triangle's signed distances against the other plane independently; parallel
offset planes are now correctly reported as non-intersecting, and its original
test still passes.

## Lean target

Not yet formalized.
