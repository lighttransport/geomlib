# Contract: 2D segment intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Orientation-based segment intersection classification (computational geometry).

## Approach

Classify endpoint orientations, detect proper crossings, endpoint touches and
collinear overlaps, and compute intersection parameters on both segments by
linear interpolation.

## API / constraints

`geomlib_segment2d_intersect(a0,a1,b0,b1)` returns a kind, crossing point or
collinear overlap, and parameters on both segments. C11, wasm32 safe.

## Acceptance

`tests/test_segment2d.cpp`: crossing segments meet at (0,0) with parameters
0.5/0.5; endpoint touching is classified; collinear overlap [1,2] returned;
parallel disjoint gives NONE. WASM smoke checks a crossing.

## Lean target

Not yet formalized.
