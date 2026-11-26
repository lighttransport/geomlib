# Contract: 2D circle–AABB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard closest-point test between a circle and an axis-aligned rectangle.

## Approach

Clamp the circle centre to the rectangle, compute the Euclidean distance to the
closest point and subtract the radius; negative values indicate penetration.

## API / constraints

`geomlib_circle_aabb_overlap(center, radius, bmin, bmax, out)` requires a
non-negative radius. C11, wasm32 safe.

## Acceptance

`tests/test_circle_aabb.cpp`: a centre one unit from the face misses by 0.5 at
radius 0.5 and penetrates by 0.5 at radius 1.5; an interior centre penetrates
by the radius. WASM smoke checks the penetrating case.

## Lean target

Not yet formalized.
