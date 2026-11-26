# Contract: AABB–AABB distance

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard per-axis gap decomposition for axis-aligned boxes.

## Approach

On each axis compute the gap between the intervals (zero when they overlap),
choose the closest coordinates on both boxes, and combine the gaps by
Pythagoras. Overlapping boxes yield distance zero and a shared point in the
intersection.

## API / constraints

`geomlib_aabb_distance(amin, amax, bmin, bmax, out)` returns 1 and fills the
distance and closest points. C11, wasm32 safe.

## Acceptance

`tests/test_aabb_distance.cpp`: unit boxes one apart on x give distance 1 with
closest points (1,0.5,0.5) and (2,0.5,0.5); overlapping boxes give 0; a
diagonal offset gives √2. WASM smoke checks the unit gap.

## Lean target

Not yet formalized.
