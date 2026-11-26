# Contract: 2D point-in-polygon

**Status:** implemented (Wave 3).

## Source (paper only)

- Ray-crossing parity and winding-number algorithms for simple polygons.

## Approach

`geomlib_point_in_polygon` uses ray crossing; `geomlib_polygon_winding`
computes the integer winding number by counting upward/downward edge crossings;
`geomlib_polygon_signed_distance` returns minus the distance to the closest
edge for interior points and positive otherwise.

## API / constraints

All take a simple polygon and a query point. Points exactly on the boundary
are treated inconsistently by ray crossing; use the signed distance for robust
classification near edges. C11, wasm32 safe.

## Acceptance

`tests/test_point_in_polygon.cpp`: a unit square contains (0.2, 0.3) with
winding +1 and signed distance −0.7; (2, 2) is outside with distance √2; a CW
square contains the point with winding −1. WASM smoke checks both cases.

## Lean target

Not yet formalized.
