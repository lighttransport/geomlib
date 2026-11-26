# Contract: 2D circle–circle intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic radical-line intersection of two circles in the plane.

## Approach

Compute the centre distance, classify none / tangent / two points / contained /
coincident, then solve `a = (r0² − r1² + d²) / 2d` and
`h = sqrt(r0² − a²)` to obtain the two symmetric intersection points.

## API / constraints

`geomlib_circle2d_intersect(c0, r0, c1, r1)` returns a kind and up to two points.
C11, wasm32 safe.

## Acceptance

`tests/test_circle_circle.cpp`: unit circles one apart intersect at
(0.5, ±√3/2); circles two apart are tangent at (1,0); three apart are disjoint;
containment and coincidence are classified. WASM smoke checks the two-point case.

## Lean target

Not yet formalized.
