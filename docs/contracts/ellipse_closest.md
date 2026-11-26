# Contract: closest point on an ellipse

**Status:** implemented (Wave 3).

## Source (paper only)

- Newton solution of the stationary-distance condition on a parametric ellipse.

## Approach

For parameter t, the ellipse point is `center + R(angle)*(a cos t, b sin t)`.
Solve `(E(t) - q) . E'(t) = 0` with Newton iteration, using eight starting
parameters and a clamped step to avoid local minima, then keep the best point.

## API / constraints

`geomlib_ellipse_closest_point(q, center, a, b, angle, out)` returns the closest
boundary point, distance and parameter. C11, wasm32 safe.

## Acceptance

`tests/test_ellipse_closest.cpp`: a query at (3,0) on semi-axes (2,1) projects
to (2,0) at distance 1; (0,2) projects to (0,1); the centre projects to the
nearer minor-axis point; a rotated major-axis query is exact. WASM smoke checks
the first case.

## Lean target

Not yet formalized.
