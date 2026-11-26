# Contract: ellipse point queries

**Status:** implemented (Wave 3).

## Source (paper only)

- Implicit equation of an oriented ellipse in the plane.

## Approach

Rotate the query into the ellipse's principal frame by −angle and evaluate
`(x/a)^2 + (y/b)^2`, which is at most 1 inside or on the boundary. Uses the
same ellipse representation as `geomlib_fit_ellipse`.

## API / constraints

`geomlib_ellipse_implicit(q, center, a, b, angle)` returns the normalized value;
`geomlib_point_in_ellipse` returns the containment flag. C11, wasm32 safe.

## Acceptance

`tests/test_ellipse.cpp`: an ellipse with semi-axes (2,1) gives implicit value
0.25 at (1,0); the minor-axis point (0,2) is outside while the major-axis
boundary point (2,0) is inside; rotated points behave accordingly. WASM smoke
checks the implicit value.

## Lean target

Not yet formalized.
