# Contract: iterative Laplacian curve smoothing

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard discrete Laplacian smoothing of polylines.

## Approach

Jacobi iteration replacing each interior vertex with a blend toward the average
of its two neighbours, with both endpoints pinned.

## API / constraints

`geomlib_smooth_polyline(points, n, iterations, lambda)` modifies points in
place; `lambda` is clamped to [0,1]. C11, wasm32 safe.

## Acceptance

`tests/test_curve_smooth.cpp`: a 21-point zigzag's roughness drops from 3.04 to
0.000029 after 100 iterations with endpoints fixed; a straight line is
unchanged. WASM smoke checks reduced roughness.

## Lean target

Not yet formalized.
