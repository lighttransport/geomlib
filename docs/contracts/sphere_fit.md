# Contract: 3D least-squares sphere fitting

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard algebraic least-squares sphere fit (linearized sphere equation).

## Approach

Solve the linear system for coefficients (a, b, c, d) in
`x²+y²+z² + a x + b y + c z + d = 0` by Gaussian elimination on the normal
equations, then recover the centre `(-a/2, -b/2, -c/2)` and radius from the
constant term.

## API / constraints

`geomlib_fit_sphere(points, n, out)` requires at least 4 points and a
non-degenerate configuration. C11, wasm32 safe.

## Acceptance

`tests/test_sphere_fit.cpp`: 400 points on a radius-2 sphere centred at
(1, -2, 3) recover the centre and radius with zero residual. WASM smoke checks
a unit sphere.

## Lean target

Not yet formalized.
