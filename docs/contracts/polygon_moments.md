# Contract: 2D polygon moments

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard shoelace formulas for polygon area, centroid and second moments.

## Approach

Accumulate cross products over polygon edges to obtain signed area, first
moments and second moments about the origin, then shift the second moments to
the centroid with the parallel-axis theorem. Winding direction is normalized.

## API / constraints

`geomlib_compute_polygon_moments(polygon, n, out)` for a simple polygon. C11,
wasm32 safe.

## Acceptance

`tests/test_polygon_moments.cpp`: a unit square gives area 1.0000, centroid
(0.5, 0.5) and inertia (1/12, 0, 1/12); a 2×1 rectangle gives 2.0000 and
(1/6, 2/3). Reversed winding leaves the results unchanged. WASM smoke checks
the square.

## Lean target

Not yet formalized.
