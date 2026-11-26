# Contract: 2D constrained Delaunay triangulation

**Status:** implemented (Wave 3).

## Source (paper only)

- Chew, "Constrained Delaunay Triangulations", Algorithmica 1989.
- De Floriani, Falcidieno, Nagy, Pienovi, "Delaunay-based representation of
  surfaces defined over arbitrarily shaped domains", 1989.

## Approach

Start from Bowyer-Watson Delaunay. For each constraint segment, repeatedly find
a crossed edge and flip it if the quadrilateral is convex; otherwise insert a
Steiner point at the crossing and split the two incident triangles. Then restore
the empty-circumcircle property by flipping all non-constrained violating edges.

## API / constraints

`geomlib_cdt2d_build(points, n, segments, m, out)`,
`geomlib_cdt2d_has_edge`, `geomlib_cdt2d_free`. Crossing constraints are allowed
but resolved by Steiner insertion, so the constraint endpoints may be joined by
a chain rather than a single edge. C11, wasm32 safe.

## Acceptance

`tests/test_cdt2d.cpp`: a square with diagonal (0,2) gives 2 triangles, area 1,
and the diagonal present; a constraint crossing interior Delaunay edges is
recovered with area preserved; all four square rim edges are present. WASM smoke
checks the diagonal case.

## Lean target

Not yet formalized.
