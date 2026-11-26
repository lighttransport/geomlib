# Contract: convex polygon minimum width

**Status:** implemented (Wave 3).

## Source (paper only)

- Shamos, "Computational Geometry", 1978; rotating calipers for antipodal
  vertex-edge pairs.

## Approach

For each edge, advance the antipodal vertex while its perpendicular distance to
the edge increases, and take the minimum of those distances. Runs in O(n) for a
convex polygon in consistent order.

## API / constraints

`geomlib_convex_polygon_min_width(poly, n)` returns 0 for fewer than three
vertices. C11, wasm32 safe.

## Acceptance

`tests/test_polygon_width.cpp`: a unit square gives 1, a 4×2 rectangle gives 2,
an equilateral triangle of side 2 gives √3, and a 128-gon radius 2 gives
3.998795. WASM smoke checks the unit square.

## Lean target

Not yet formalized.
