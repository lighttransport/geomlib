# Contract: circle–ellipse overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-point collision test between a circle and an oriented ellipse.

## Approach

Project the circle centre onto the ellipse boundary with the existing closest-
point solver and subtract the radius. Centre containment is reported separately
because a small circle centred inside the ellipse overlaps even when its
boundary gap is positive.

## API / constraints

`geomlib_circle_ellipse_overlap(circle_center, radius, ellipse_center, a, b,
angle, out)` uses the existing ellipse representation. C11, wasm32 safe.

## Acceptance

`tests/test_circle_ellipse.cpp`: a circle one unit from semi-axes (2,1) misses
with radius 0.5 and penetrates by 0.5 with radius 1.5; a centred circle overlaps
via containment. WASM smoke checks the penetrating case.

## Lean target

Not yet formalized.
