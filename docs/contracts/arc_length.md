# Contract: polyline arc length

**Status:** implemented (Wave 3).

## Source (paper only)

- Cumulative chord-length parameterization of a polyline.

## Approach

Build a cumulative distance table in double precision, provide the total length,
and linearly interpolate the polyline at a requested arc length with clamping
at both ends.

## API / constraints

`geomlib_polyline_arc_table`, `geomlib_polyline_length` and
`geomlib_polyline_at_length` operate on 3D polylines. C11, wasm32 safe.

## Acceptance

`tests/test_arc_length.cpp`: an L-shaped polyline has cumulative lengths
(0,3,7), total 7; queries at s=1.5 and s=5 give (1.5,0,0) and (3,2,0), and
over-range queries clamp to the end. WASM smoke checks the total length.

## Lean target

Not yet formalized.
