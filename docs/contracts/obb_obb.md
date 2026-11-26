# Contract: OBB–OBB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Gottschalk, "Separating Axis Theorem for Oriented Bounding Boxes", 2000.

## Approach

Test separation over 15 axes: the 3 axes of each box and the 9 cross products.
For each axis compare the centre separation with the sum of projected radii
from both boxes. Degenerate axes are skipped.

## API / constraints

`geomlib_obb_obb_overlap(a, b)` returns 1 on overlap. C11, wasm32 safe.

## Acceptance

`tests/test_obb_obb.cpp`: coincident OBBs overlap, a translated box is separated,
and a 45-degree rotated overlapping box is detected. WASM smoke checks the
separated case.

## Lean target

Not yet formalized.
