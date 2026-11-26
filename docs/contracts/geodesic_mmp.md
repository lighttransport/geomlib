# Contract: continuous geodesic via MMP window propagation

**Status:** implemented (Wave 3).

## Source (paper only)

- Mitchell, Mount and Papadimitriou, "The Continuous Dijkstra Problem", SIAM
  Journal on Computing 1987.

## Approach

Continuous Dijkstra by window propagation. A window stores an interval on a
mesh edge plus the source unfolded into the canonical plane of a face incident
to that edge. Windows are popped in increasing minimum distance; crossing the
shared edge mirrors the source, and the reachable interval is propagated to the
other two edges by intersecting the source fan with them. This yields a true
continuous surface distance rather than the face-strip upper bound of
`geomlib_geodesic_length_unfolding`.

## API / constraints

`geomlib_geodesic_mmp(mesh, source, target, length)` works on triangle meshes
with valid vertex indices. A safety cap limits window chains to 64 face
crossings and 40,000 windows, so extremely long wrapping paths may fall back to
the best found so far. C11, wasm32 safe.

## Acceptance

`tests/test_geodesic_mmp.cpp`: a flat zigzag strip from (0,0) to (2,1) gives
sqrt(5) even though no face contains the vertex triplet of the edge shortest
path; opposite corners of the unit cube give sqrt(5). WASM smoke checks the
cube case.

## Lean target

Not yet formalized.
