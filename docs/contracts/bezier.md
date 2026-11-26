# Contract: Bezier curve evaluation and subdivision

**Status:** implemented (Wave 3).

## Source (paper only)

- de Casteljau, "Mathematische Algorithmen: Eine Einführung", 1963; the
  recursive linear-interpolation construction of Bezier curves.

## Approach

Build the full de Casteljau interpolation triangle in double precision. The
single remaining point is the curve value; the first entry of each level forms
the left subcurve and the last entry of each reversed level forms the right
subcurve. The derivative is a Bezier curve over scaled control-point differences.

## API / constraints

`geomlib_bezier_evaluate`, `geomlib_bezier_derivative` and
`geomlib_bezier_subdivide` accept `n` control points (degree n-1) and clamp t to
[0,1]. C11, wasm32 safe.

## Acceptance

`tests/test_bezier.cpp`: endpoints are interpolated; linear controls give
x=1.5 at t=0.5; derivative at t=0 is 3(P1-P0); subdivision halves join at the
evaluated point and reproduce the original curve. WASM smoke checks the midpoint.

## Lean target

Not yet formalized.
