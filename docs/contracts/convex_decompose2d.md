# Contract: 2D convex polygon decomposition

**Status:** implemented (Wave 3).

## Source (paper only)

- Chazelle and Dobkin, "Decomposing a Polygon into Its Convex Parts", 1979.

## Approach

If the polygon is convex, return it. Otherwise score every valid diagonal from
a reflex vertex by the total number of reflex vertices in the two resulting
sub-polygons, keep the best cut, and recurse on both parts up to a depth limit.

## API / constraints

`geomlib_convex_decompose2d(polygon, n, out)` and
`geomlib_poly_decomposition_free`. Input must be a simple polygon. The result
is a greedy decomposition, not necessarily minimum-cardinality. C11, wasm32
safe.

## Acceptance

`tests/test_convex_decompose2d.cpp`: an L-shape splits into two convex parts
with total area 3.0000; a square stays a single part. WASM smoke checks two
parts for the L-shape.

## Lean target

Not yet formalized.
