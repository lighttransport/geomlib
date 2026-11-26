# Contract: convex polygon diameter

**Status:** implemented (Wave 3).

## Source (paper only)

- Shamos, "Computational Geometry", 1978; rotating calipers for antipodal pairs.

## Approach

Advance an antipodal vertex j while the triangle area formed with each edge
increases, evaluating the distance at every visited pair. Runs in O(n) for a
convex polygon in consistent order.

## API / constraints

`geomlib_convex_polygon_diameter(poly, n)` returns the maximum vertex distance;
0 for fewer than two vertices. C11, wasm32 safe.

## Acceptance

`tests/test_polygon_diameter.cpp`: a unit square gives √2, a 4×2 rectangle gives
√20, and a 128-gon radius 2 gives 4.000000. WASM smoke checks the square.

## Lean target

Not yet formalized.
