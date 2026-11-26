# Contract: convex polygon clipping

**Status:** implemented (Wave 3).

## Source (paper only)

- Sutherland and Hodgman, "Reentrant Polygon Clipping", CACM 1974.

## Approach

Clip the subject polygon sequentially against each edge half-plane of the
convex clip polygon, inserting intersection points where the inside/outside
status changes along an edge.

## API / constraints

`geomlib_clip_polygon_convex(subject, ns, clip, nc, out)` fills an existing
`geomlib_polygon2d`; an empty result means no intersection. Both inputs should
be CCW and the clip polygon convex. C11, wasm32 safe.

## Acceptance

`tests/test_polygon_clip.cpp`: clipping a square by itself gives 4 vertices
and area 4.0000; a half-plane gives area 2.0000; a disjoint box gives an empty
polygon. WASM smoke checks the half-plane case.

## Lean target

Not yet formalized.
