# Contract: convex polygon intersection area

**Status:** implemented (Wave 3).

## Source (paper only)

- Sutherland–Hodgman polygon clipping and the shoelace area formula.

## Approach

Clip the first convex polygon against the second with the existing clipper,
then compute the absolute signed area of the resulting polygon. Returns 0 for
disjoint polygons.

## API / constraints

`geomlib_polygon_intersection_area(a, na, b, nb)` assumes both polygons are
convex and consistently wound. C11, wasm32 safe.

## Acceptance

`tests/test_polygon_intersection_area.cpp`: half-overlapping unit squares give
0.25; identical squares give 1; disjoint squares give 0. WASM smoke checks the
half-overlap case.

## Lean target

Not yet formalized.
