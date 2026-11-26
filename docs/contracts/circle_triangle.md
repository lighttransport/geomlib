# Contract: 2D circle–triangle intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard circle-triangle overlap test via containment and closest edge point.

## Approach

Test centre containment with the polygon point test, vertex containment with
the radius, and compute the closest boundary point from the three edges. The
signed distance is the centre-to-boundary distance minus the radius; overlap
holds if the centre is inside, a vertex is inside, or the signed distance is
non-positive.

## API / constraints

`geomlib_circle_triangle_intersect(center, radius, a, b, c, out)` requires a
non-negative radius. C11, wasm32 safe.

## Acceptance

`tests/test_circle_triangle.cpp`: an interior centre intersects even with a
positive edge gap; a distant circle misses; a large circle contains the
triangle; a circle tangent to an edge gives signed distance zero at (1,0).
WASM smoke checks interior containment.

## Lean target

Not yet formalized.
