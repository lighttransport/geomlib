# Contract: capsule–AABB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard segment-to-box closest-point enumeration for capsule collision.

## Approach

Test whether the axis segment crosses the box by a slab test. Otherwise take
the minimum of the endpoint-to-box distances and the 12 box-edge-to-segment
closest distances. Subtract the capsule radius for the signed penetration depth.

## API / constraints

`geomlib_segment_aabb_distance(a,b,bmin,bmax)` exposes the raw axis distance;
`geomlib_capsule_aabb_overlap(a,b,r,bmin,bmax,out)` returns the overlap flag.
C11, wasm32 safe.

## Acceptance

`tests/test_capsule_aabb.cpp`: a crossing axis gives penetration 0.2; a segment
one unit from the box gives distance 1; radii 0.5 and 1.5 respectively miss and
penetrate by 0.5. WASM smoke checks the penetrating case.

## Lean target

Not yet formalized.
