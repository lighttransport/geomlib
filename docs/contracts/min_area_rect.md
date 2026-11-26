# Contract: minimum-area bounding rectangle

**Status:** implemented (Wave 3).

## Source (paper only)

- Freeman and Shapira, "Finding the Minimum-Area Rectangular Enclosure of a
  Polygon", Pattern Recognition 1975 (rotating calipers).

## Approach

Compute the convex hull, then for each hull edge use that edge direction and
its perpendicular as rectangle axes and measure the projected extents of all
hull vertices. Keep the orientation with minimum area.

## API / constraints

`geomlib_minimum_area_rectangle(points, n, out)` returns the centre, two unit
axes, half-extents and area. The implementation is O(h²) in hull size rather
than the optimal O(h). C11, wasm32 safe.

## Acceptance

`tests/test_min_area_rect.cpp`: a unit square gives area 1.0000 and centre
(0.5, 0.5); a 30°-rotated 2×1 rectangle gives area 2.0000; a right triangle
gives area 1.0000. All points are inside. WASM smoke checks a square.

## Lean target

Not yet formalized.
