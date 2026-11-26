# Contract: polyline simplification

**Status:** implemented (Wave 3).

## Source (paper only)

- Douglas and Peucker, "Algorithms for the Reduction of the Number of Points
  Required to Represent a Digitized Line or its Caricature", 1973.

## Approach

Keep the two endpoints, find the interior point farthest from the segment
joining them, and recurse on both halves while the distance exceeds epsilon.
`geomlib_simplify_max_error` reports the maximum distance from any input point
to the simplified polyline.

## API / constraints

`geomlib_simplify_polyline(points, n, epsilon, indices, out_count)` and
`geomlib_simplify_max_error`. The output indices are in original order. For
closed polylines, callers should pass each half separately. C11, wasm32 safe.

## Acceptance

`tests/test_simplify.cpp`: 50 collinear points reduce to 2 endpoints with zero
error; a sparse-peak zigzag reduces 21 to 14 points with error within tolerance;
a triangle keeps all three corners. WASM smoke checks a 24-point polyline is
reduced.

## Lean target

Not yet formalized.
