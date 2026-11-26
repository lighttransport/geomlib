# Contract: continuous geodesic length via unfolding

**Status:** implemented (Wave 3).

## Source (papers only)

- Mitchell, Mount, Papadimitriou, "The Continuous Dijkstra Problem", 1987
  (unfolding/straightest path concept).
- Kimmel & Sethian, "Computing Geodesic Paths on Manifolds", PNAS 1998.

## Approach

Find a shortest edge path, then unfold each traversed triangle into a plane by
solving a two-circle intersection across the shared edge. The straight-line
distance between the unfolded endpoints is the exact shortest path *within
that face strip*. For arbitrary face sequences found by graph Dijkstra this is
an upper bound on the true continuous geodesic; a full MMP sequence-tree
implementation is future work.

## API / constraints

`geomlib_geodesic_length_unfolding(mesh, source, target, length)`.
C11, wasm32 safe. Zero-length duplicate vertices are skipped.

## Acceptance

`tests/test_geodesic_exact.cpp`: unit-square diagonal unfolds to sqrt(2)
exactly; a sphere edge path reproduces its edge length. WASM smoke checks a
plane path length.

## Lean target

Not yet formalized.
