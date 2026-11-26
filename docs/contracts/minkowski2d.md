# Contract: 2D convex Minkowski sum

**Status:** implemented (Wave 3).

## Source (paper only)

- Minkowski, "Volumen und Oberfläche", 1903; the convex-polynomial identity
  `conv(A) + conv(B) = conv{a+b | a in A, b in B}` is the classical basis for
  computing sums by vertex pairs.

## Approach

Form all `na*nb` pairwise vertex sums and return their convex hull. This is
exact when both inputs are convex and yields the Minkowski sum of the convex
hulls otherwise.

## API / constraints

`geomlib_minkowski_sum2d(a, na, b, nb, out)`, returning a `geomlib_hull2d`
owned by the caller. Cost is O(na*nb log(na*nb)). C11, wasm32 safe.

## Acceptance

`tests/test_minkowski2d.cpp`: two unit squares give a 2×2 square of area 4;
a square plus a point stays area 1; triangle plus square gives a convex
pentagon of area 3.5. WASM smoke checks the square case.

## Lean target

Not yet formalized.
