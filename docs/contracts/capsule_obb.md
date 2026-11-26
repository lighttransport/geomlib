# Contract: capsule–OBB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Local-frame closest-point collision between a capsule and an oriented box.

## Approach

Transform both capsule axis endpoints into the OBB's orthonormal frame, compute
the segment-to-box distance there with the shared segment-AABB routine, and
subtract the capsule radius for signed penetration depth.

## API / constraints

`geomlib_capsule_obb_overlap(a, b, r, obb, out)` uses the existing `geomlib_obb`
type and requires a non-negative radius. C11, wasm32 safe.

## Acceptance

`tests/test_capsule_obb.cpp`: an axis one unit from a unit box misses with
radius 0.5 and penetrates by 0.5 with radius 1.5; a rotated box gives the same
gap in local coordinates. WASM smoke checks the penetrating case.

## Lean target

Not yet formalized.
