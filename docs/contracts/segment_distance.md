# Contract: 3D segment-segment closest points

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard clamped parametric closest-point solution for two line segments.

## Approach

Solve the infinite-line closest-point system, clamp the parameters to [0,1],
and handle parallel/degenerate axes explicitly. Returns both closest points,
their segment parameters and the distance.

## API / constraints

`geomlib_segment_closest_points(a0, a1, b0, b1, out)` for finite segments. C11,
wasm32 safe.

## Acceptance

`tests/test_segment_distance.cpp`: crossing segments give 0, parallel segments
separated by 2 give 2, and skew segments give 1 with parameters in [0,1]. WASM
smoke checks a known 3-unit separation.

## Lean target

Not yet formalized.
