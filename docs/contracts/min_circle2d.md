# Contract: exact 2D minimum enclosing circle

**Status:** implemented (Wave 3).

## Source (paper only)

- Welzl, "Smallest Enclosing Disks (Balls and Ellipsoids)", 1991.

## Approach

Randomized incremental construction with support sets of up to three points.
Trivial circles from all nonempty subsets of the support set are enumerated and
the smallest one containing the whole support set is kept.

## API / constraints

`geomlib_minimum_enclosing_circle(points, n, out)` for up to 512 points. C11,
wasm32 safe.

## Acceptance

`tests/test_min_circle2d.cpp`: 200 points on a unit circle give radius 1.000000;
two points 4 apart give radius 2.0000; an obtuse triangle gives the longest-edge
diameter 2.0000. WASM smoke checks the unit circle.

## Lean target

Not yet formalized.
