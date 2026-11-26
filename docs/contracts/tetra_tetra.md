# Contract: tetrahedron–tetrahedron intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-feature enumeration for two convex tetrahedra.

## Approach

Test vertex containment in both directions and all 6 tetrahedron edges against
all 4 opposing faces using the existing segment-triangle predicate. Any positive
case implies intersection.

## API / constraints

`geomlib_tetrahedra_intersect(a0,a1,a2,a3,b0,b1,b2,b3)` returns 1 on overlap.
C11, wasm32 safe.

## Acceptance

`tests/test_tetra_tetra.cpp`: identical, contained and partially overlapping
tetrahedra intersect; a translated disjoint tetrahedron does not. WASM smoke
checks partial overlap.

## Lean target

Not yet formalized.
