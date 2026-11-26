# Contract: Chaikin polyline subdivision

**Status:** implemented (Wave 3).

## Source (paper only)

- Chaikin, "An Algorithm for High Speed Curve Generation", Computer Graphics
  and Image Processing 1974.

## Approach

Repeatedly replace each segment with two points at 1/4 and 3/4 along it. For
open curves the endpoints are preserved explicitly. The output count is
`(count-1)*2^iterations + 1`.

## API / constraints

`geomlib_chaikin_subdivide(points, n, iterations, out, out_capacity)` returns
the number written or 0 if the capacity is insufficient. C11, wasm32 safe.

## Acceptance

`tests/test_chaikin.cpp`: two levels on a straight 3-point line give 9 points,
all collinear with endpoints 0 and 6; one level on a corner gives 5 points with
endpoints preserved and the apex rounded. WASM smoke checks a 5-point result.

## Lean target

Not yet formalized.
