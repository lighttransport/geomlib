# Contract: 2D line–line intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Cross-product solution of two parametric lines in the plane.

## Approach

Let `u = b-a`, `v = d-c`, `w = c-a`. If `cross(u,v)` is non-zero, the parameter
along the first line is `t = cross(w,v)/cross(u,v)` and the intersection point
is `a + t u`. A zero denominator means parallel or collinear, classified via
`cross(w,u)`.

## API / constraints

`geomlib_line2d_intersect(a,b,c,d)` uses infinite lines and returns none / point
/ collinear with the intersection parameters. C11, wasm32 safe.

## Acceptance

`tests/test_line2d.cpp`: lines (0,0)-(2,0) and (1,-1)-(1,1) intersect at (1,0)
with t=u=0.5; parallel distinct lines miss; collinear lines are classified.
WASM smoke checks the crossing case.

## Lean target

Not yet formalized.
