# Contract: 2D least-squares circle fitting

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard algebraic least-squares circle fit (linearized circle equation).

## Approach

Solve the linear system for coefficients (a, b, c) in
`x² + y² + a x + b y + c = 0` by Gaussian elimination on the normal equations,
then recover the centre `(-a/2, -b/2)` and radius from the constant term.

## API / constraints

`geomlib_fit_circle(points, n, out)` requires at least 3 points and a
non-degenerate configuration. C11, wasm32 safe.

## Acceptance

`tests/test_circle_fit.cpp`: 200 points on a circle of radius 3 centred at
(2, -1) recover the centre and radius with zero residual. WASM smoke checks a
radius-2 circle.

## Lean target

Not yet formalized.
