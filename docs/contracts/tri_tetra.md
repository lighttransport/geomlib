# Contract: triangle–tetrahedron intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-feature enumeration between a triangle and a tetrahedron.

## Approach

Return true when a triangle vertex lies in the tetrahedron, a tetrahedron
vertex lies in the triangle, a triangle edge meets a tetrahedron face, or a
tetrahedron edge meets the triangle. Reuses the existing point-in-tetrahedron
and segment-triangle helpers.

## API / constraints

`geomlib_triangle_tetrahedron_intersect(a0,a1,a2,v0,v1,v2,v3)` returns 1 on
intersection. C11, wasm32 safe.

## Acceptance

`tests/test_tri_tetra.cpp`: a cutting triangle, an interior triangle and a
shared face all intersect; a distant triangle misses. WASM smoke checks the
cutting case.

## Note

`geomlib_segment_triangle_intersect` was corrected to work as a pure predicate
with a NULL output pointer; its original test still passes.

## Lean target

Not yet formalized.
