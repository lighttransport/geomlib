# Contract: Catmull-Rom spline resampling

**Status:** implemented (Wave 3).

## Source (paper only)

- Catmull and Rom, "A Class of Local Interpolating Splines", 1974.
- Barry and Goldman, "A Recursive Evaluation Algorithm for a Class of Catmull-
  Rom Curves", SIGGRAPH 1988 (pyramid evaluation).

## Approach

For each control-point segment, build knots from the chord lengths (square
root for centripetal, uniform spacing otherwise) and evaluate the Barry-Goldman
pyramid. Coincident controls get epsilon-separated knots. Closed curves repeat
the first point at the end.

## API / constraints

`geomlib_resample_catmull_rom(points, n, closed, samples_per_segment, type, out)`
and `geomlib_catmull_output_count`. The caller allocates `out`. C11, wasm32 safe.

## Acceptance

`tests/test_spline.cpp`: open collinear controls give exactly 25 samples on the
line with interpolated endpoints and control points; a closed diamond gives a
loop with radii in [0.884, 1.000]; uniform and centripetal differ by 2.23 on
unevenly spaced controls. WASM smoke checks a 25-point centripetal curve.

## Lean target

Not yet formalized.
