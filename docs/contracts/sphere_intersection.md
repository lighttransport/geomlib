# Contract: sphere–sphere intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic intersection of two spheres in 3D Euclidean space.

## Approach

Compute the axis from the centre difference, find the circle plane position by
the radical-plane formula, and derive the radius from Pythagoras. Classify no
intersection, single-point tangency, a circle, or coincident spheres.

## API / constraints

`geomlib_sphere_sphere_intersection(c0, r0, c1, r1)` returns kind, circle
centre, unit axis and radius. C11, wasm32 safe.

## Acceptance

`tests/test_sphere_intersection.cpp`: two radius-2 spheres 3 apart give a
circle at x=1.5 with radius √1.75 ≈ 1.3229; external tangency, disjoint and
coincident cases classified. WASM smoke checks the circle radius.

## Lean target

Not yet formalized.
