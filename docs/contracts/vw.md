# Contract: Visvalingam-Whyatt polyline simplification

**Status:** implemented (Wave 3).

## Source (paper only)

- Visvalingam and Whyatt, "Line Generalisation by Repeated Elimination of
  Points", Cartographic Information Systems 1993.

## Approach

The effective area of an interior point is the area of the triangle with its
current neighbours. Repeatedly remove the smallest-area point while that area
is at most `area_epsilon`, updating the two neighbours' areas after each
removal. Endpoints are never removed.

## API / constraints

`geomlib_simplify_polyline_vw(points, n, area_epsilon, indices, out_count)`
returns kept indices in original order. `geomlib_vw_effective_area` exposes the
per-point area. O(n²) worst case. C11, wasm32 safe.

## Acceptance

`tests/test_vw.cpp`: an isolated spike has effective area 5 exactly; simplifying
with epsilon 0.1 keeps endpoints, the spike and its ramp points (5 kept); a
straight line keeps only its endpoints. WASM smoke checks the 5-kept result.

## Lean target

Not yet formalized.
