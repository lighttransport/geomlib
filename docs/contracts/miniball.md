# Contract: minimum enclosing sphere

**Status:** implemented (Wave 3).

## Source (paper only)

- Welzl, "Smallest Enclosing Disks (Balls and Ellipsoids)", New Results and
  New Trends in Computer Graphics, 1991.

## Approach

Randomized incremental construction: recursively solve the subproblem without
a point; if that point lies outside the result, add it to the support set.
Trivial balls from up to four support points are found exactly by enumerating
subsets and solving the corresponding linear system for the centre.

## API / constraints

`geomlib_minimum_enclosing_sphere(points, n, out)` for up to 512 points. Uses
the shared `geomlib_bounding_sphere` result type. C11, wasm32 safe.

## Acceptance

`tests/test_miniball.cpp`: 200 points on a unit sphere give radius 1.00000
(Ritter gives 1.03182); a unit box's corners give exactly 1.73205. WASM smoke
checks the box radius.

## Lean target

Not yet formalized.
