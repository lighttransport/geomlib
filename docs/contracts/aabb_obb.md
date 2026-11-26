# Contract: AABB–OBB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Separating-axis test for two oriented boxes, with the AABB treated as a box
  whose axes are the world basis.

## Approach

Test separation over 15 axes: the 3 world axes, the 3 OBB axes and the 9 cross
products. For each axis compare the centre separation with the sum of projected
box radii. Degenerate axes are skipped.

## API / constraints

`geomlib_aabb_obb_overlap(amin, amax, obb)` returns 1 on overlap. C11, wasm32
safe.

## Acceptance

`tests/test_aabb_obb.cpp`: coincident AABB/OBB overlap, a translated box is
separated, and a rotated overlapping box is detected. WASM smoke checks the
separated case.

## Lean target

Not yet formalized.
