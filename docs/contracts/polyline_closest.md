# Contract: closest point on a polyline

**Status:** implemented (Wave 3).

## Source (paper only)

- Finite-segment projection with endpoint clamping.

## Approach

Project the query onto each segment, clamp the parameter to [0,1], and keep
the nearest projection together with the segment index and parameter.

## API / constraints

`geomlib_polyline_closest_point(points, n, q, out)` supports any `n >= 1` and
returns the closest point, distance, segment index and local parameter. C11,
wasm32 safe.

## Acceptance

`tests/test_polyline_closest.cpp`: query (1,1,0) on an L-shaped polyline gives
(1,0,0) at distance 1 on segment 0 with t=0.5; query (3,2,0) clamps to the end
(2,2,0) on segment 1 with t=1; a single-point polyline gives distance 5. WASM
smoke checks the first case.

## Lean target

Not yet formalized.
