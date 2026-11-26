# Contract: minimum spanning tree

**Status:** implemented (Wave 3).

## Source (paper only)

- Prim, "Shortest Connection Networks and Some Generalizations", Bell System
  Technical Journal 1957.

## Approach

Grow a tree from the first point, repeatedly adding the nearest unvisited point
by squared Euclidean distance, accumulating actual edge lengths.

## API / constraints

`geomlib_compute_mst(points, n, dim, out)` for 2D or 3D points. O(n²) time.
C11, wasm32 safe.

## Acceptance

`tests/test_mst.cpp`: a unit square gives 3 edges, total weight 3.0000 and a
connected tree; collinear points 0, 1, 3 give weight 3. WASM smoke checks the
square.

## Lean target

Not yet formalized.
