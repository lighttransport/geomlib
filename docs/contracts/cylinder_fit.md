# Contract: least-squares cylinder fitting

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard nonlinear least-squares fit of an infinite cylinder to 3D points.

## Approach

Initialize the axis from a caller-provided direction and the axis point from
the centroid, then minimize squared radial residuals by finite-difference
gradient descent with backtracking. The axis is renormalized after each step.

## API / constraints

`geomlib_fit_cylinder(points, n, initial_axis, out)` requires at least 5 points
and a reasonable initial axis direction. C11, wasm32 safe.

## Acceptance

`tests/test_cylinder_fit.cpp`: 300 points on a radius-1.5 cylinder centred at
(2, −1) along z recover the axis, centre, radius and zero residual. WASM smoke
checks a radius-2 cylinder.

## Lean target

Not yet formalized.
