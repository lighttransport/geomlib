# Contract: triangle-triangle intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Möller, "A Fast Triangle-Triangle Intersection Test", JGT 1997; this
  implementation uses bidirectional edge/triangle tests with a coplanar
  2D edge-pair fallback.

## Approach

Trace each edge of one triangle against the other triangle, collect intersection
points from both directions, deduplicate, and return the most distant pair as
the intersection segment. Coplanar triangles are handled by projecting edge
pairs onto a 2D basis and solving segment intersections directly.

## API / constraints

`geomlib_triangle_intersection(a0,a1,a2,b0,b1,b2)` returns an `intersects` flag
and a segment. Coplanar overlapping faces may return a representative segment.
C11, wasm32 safe.

## Acceptance

`tests/test_tri_tri.cpp`: crossing coplanar triangles give segment
(0,−1)–(0,1); disjoint triangles do not intersect; perpendicular triangles
intersect. WASM smoke checks the coplanar case.

## Lean target

Not yet formalized.
