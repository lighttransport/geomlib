# Contract: triangle–OBB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Local-frame separating-axis test between a triangle and an oriented box.

## Approach

Transform the triangle vertices into the OBB's orthonormal frame, then reuse
the 13-axis triangle-AABB overlap test with box half-extents from the OBB.

## API / constraints

`geomlib_triangle_obb_overlap(a0,a1,a2,obb)` returns 1 on overlap. C11, wasm32
safe.

## Acceptance

`tests/test_tri_obb.cpp`: a triangle crossing a unit box overlaps, one beside it
is separated, and a rotated box is handled in local coordinates. WASM smoke
checks the crossing case.

## Lean target

Not yet formalized.
