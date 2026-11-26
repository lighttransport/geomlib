# Contract: segment–triangle intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Möller and Trumbore, "Fast, Minimum Storage Ray/Triangle Intersection",
  Journal of Graphics Tools 1997; restricted here to a finite segment.

## Approach

Apply the Möller–Trumbore solve with direction `b − a`, then accept the hit only
when the segment parameter and both barycentric coordinates lie in [0,1] within
tolerance. Coplanar segments are not reported because their intersection is not
a single point.

## API / constraints

`geomlib_segment_triangle_intersect(a, b, v0, v1, v2, out)` returns 1 on a hit
and fills t, point and barycentrics. C11, wasm32 safe.

## Acceptance

`tests/test_segment_triangle.cpp`: a segment crossing the unit right triangle at
(0.2,0.2) hits at t=0.5 with barycentrics (0.6,0.2,0.2); short, outside and
coplanar segments miss. WASM smoke checks the exact hit.

## Lean target

Not yet formalized.
