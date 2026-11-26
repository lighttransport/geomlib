# Contract: 2D Delaunay triangulation

**Status:** implemented (Wave 3).

## Source (paper only)

- Bowyer, "Computing Dirichlet Tessellations", 1981; Watson, "Computing the
  n-dimensional Delaunay tessellation", 1981.

## Approach

Incremental Bowyer-Watson insertion with a super-triangle; bad triangles whose
circumcircle contains the new point are removed, their boundary edges are
extracted, and new triangles are connected to the point. Super-triangle
triangles are dropped at the end. The exact in-circle predicate from the
Wave 0 predicates module is used in verification.

## API / constraints

`geomlib_delaunay2d(points, count, out)` returning triangle indices.
C11, wasm32 safe.

## Acceptance

`tests/test_delaunay2d.cpp`: square gives exactly 2 triangles; 40 random points
give 69 triangles satisfying the empty circumcircle property under the exact
in-circle predicate. WASM smoke triangulates a square.

## Lean target

2D Delaunay flip termination is future work (related formal work exists for
in-circle axioms in `formal/GeomLibFormal/Delaunay.lean`).
