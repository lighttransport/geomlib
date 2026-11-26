# Contract: 2D convex hull

**Status:** implemented (Wave 3).

## Source (paper only)

- Andrew, "Another Efficient Algorithm for Convex Hulls in Two Dimensions",
  Information Processing Letters 1979.

## Approach

Sort points lexicographically, build the lower chain left-to-right and the
upper chain right-to-left, popping while the last turn is non-left. Duplicate
input points are removed and collinear points are dropped from the boundary.

## API / constraints

`geomlib_convex_hull2d(points, n, out)` and `geomlib_hull2d_free`. The output
is CCW with no repeated closing point; one or two points are returned as-is.
C11, wasm32 safe.

## Acceptance

`tests/test_hull2d.cpp`: a square plus interior points gives four CCW hull
vertices, area 4.0000, and every input point inside; four collinear points give
two endpoints. WASM smoke checks the square hull.

## Lean target

Not yet formalized.
