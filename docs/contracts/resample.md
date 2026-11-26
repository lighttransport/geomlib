# Contract: arc-length polyline resampling

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard cumulative-chord-length parameterization of polylines.

## Approach

Compute cumulative segment lengths, then interpolate each requested sample at
the corresponding fraction of the total length, preserving both endpoints.

## API / constraints

`geomlib_resample_polyline(points, n, samples, out)` writes `samples` points.
Equal arc-length spacing does not imply equal Euclidean spacing across corners.
C11, wasm32 safe.

## Acceptance

`tests/test_resample.cpp`: a straight polyline with uneven original spacing
resamples to 21 points with total length 9.0000 and uniform 0.4500 spacing.
WASM smoke checks length preservation.

## Lean target

Not yet formalized.
