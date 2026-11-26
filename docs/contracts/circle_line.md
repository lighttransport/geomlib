# Contract: 2D circle–line intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic quadratic intersection of a circle with a parametric line in 2D.

## Approach

Substitute the parametric line into the circle equation, solve the quadratic,
and keep roots in [0,1] for the finite-segment variant. Coincident roots are
reported as tangency.

## API / constraints

`geomlib_circle_line_intersect` uses the infinite line through two points;
`geomlib_circle_segment_intersect` restricts to the finite segment. Both return
none / tangent / two points. C11, wasm32 safe.

## Acceptance

`tests/test_circle_line.cpp`: a secant through the unit circle gives ±1; a
tangent line touches at (0,1); a distant line misses; a segment from x=0.5 to 2
keeps only x=1. WASM smoke checks the secant.

## Lean target

Not yet formalized.
