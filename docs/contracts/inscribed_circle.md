# Contract: largest inscribed circle

**Status:** implemented (Wave 3).

## Source (paper only)

- The Chebyshev centre of a polygon is a vertex of the Voronoi diagram of its
  boundary; standard computational-geometry result.

## Approach

Sample the polygon boundary uniformly, build the clipped Voronoi diagram of the
samples, and keep the interior Voronoi vertex furthest from all polygon edges.
The result is an approximation whose accuracy depends on the sample spacing.

## API / constraints

`geomlib_largest_inscribed_circle(polygon, n, spacing, out)`. `spacing <= 0`
selects 0.02 × bounding-box diagonal. Input must be a simple polygon. C11,
wasm32 safe.

## Acceptance

`tests/test_inscribed_circle.cpp`: a square gives centre (−0.000, 0.000) and
radius 1.0000 touching at least two edges; a 4×2 rectangle gives radius 1.0000.
WASM smoke checks the square radius.

## Lean target

Not yet formalized.
