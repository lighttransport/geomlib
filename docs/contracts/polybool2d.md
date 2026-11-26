# Contract: 2D polygon Boolean operations

**Status:** implemented (Wave 3).

## Source (paper only)

- Vatti, "A Generic Solution to Polygon Clipping", SIGGRAPH 1992.
- Standard edge-arrangement Boolean construction.

## Approach

Normalize both polygons to CCW, split every edge at crossings with the other
polygon (including collinear endpoint contacts), classify each sub-edge midpoint
as inside/outside the other polygon, keep the parts selected by the Boolean
predicate, and trace the directed sub-edges into closed loops. Difference
reverses the kept edges of B so holes wind clockwise.

## API / constraints

`geomlib_polygon_boolean(a, na, b, nb, op, out)` for union, intersection and
A-minus-B, with `geomlib_poly_bool_result_free`. Inputs must be simple polygons;
coplanar overlapping edges are not merged into shared runs. C11, wasm32 safe.

## Acceptance

`tests/test_polybool2d.cpp`: overlapping 2×2 squares offset by (1,1) give areas
7 (union), 1 (intersection) and 3 (difference); a contained square gives 9, 1
and 0. WASM smoke checks the overlap area is 1.

## Lean target

Not yet formalized.
