# Contract: discrete Fréchet distance

**Status:** implemented (Wave 3).

## Source (paper only)

- Alt and Godau, "Computing the Fréchet Distance Between Two Polygonal Curves",
  IJCGA 1995; the discrete variant uses dynamic programming on the coupling.

## Approach

Fill an na×nb DP table where each entry is the minimum over the three predecessor
moves of the maximum of the previous coupling distance and the current endpoint
distance. Returns the final DP value using two rolling rows.

## API / constraints

`geomlib_frechet_distance(a, na, b, nb)` for 3D polylines; O(na·nb) time and
O(nb) memory. C11, wasm32 safe.

## Acceptance

`tests/test_frechet.cpp`: identical curves give 0; a translation by 3 gives
3.0000; the distance is symmetric. WASM smoke checks a translation by 2.

## Lean target

Not yet formalized.
