# Contract: capsule–triangle overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-feature distance between a finite segment and a triangle.

## Approach

First use the segment-triangle intersection predicate; if it hits, the distance
is zero. Otherwise enumerate endpoint-to-triangle projections and the three
triangle edges against the capsule axis, keeping the minimum. Subtract the
capsule radius for signed penetration.

## API / constraints

`geomlib_segment_triangle_distance` exposes the raw distance with witnesses;
`geomlib_capsule_triangle_overlap(a,b,r,v0,v1,v2,out)` returns the overlap flag.
C11, wasm32 safe.

## Acceptance

`tests/test_capsule_triangle.cpp`: a crossing axis penetrates by 0.2; a parallel
axis above the triangle misses with gap 1.5; an axis crossing an edge overlaps.
WASM smoke checks the penetration case.

## Lean target

Not yet formalized.
