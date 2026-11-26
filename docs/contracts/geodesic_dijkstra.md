# Contract: Dijkstra geodesic distances

**Status:** implemented (Wave 3).

## Source (paper only)

- Dijkstra, "A Note on Two Problems in Connexion with Graphs", 1959.

## Approach

Build the undirected triangle-edge graph with Euclidean edge lengths, then run
a binary-heap multi-source Dijkstra. Distances are stored in the shared
`geomlib_geodesic_distances` structure, with source indices recorded.

## API / constraints

`geomlib_geodesic_dijkstra(mesh, sources, n, out)` and
`geomlib_geodesic_dijkstra_from_vertex(mesh, source, out)`. The mesh must have
shared vertices across faces — per-face-vertex boxes are disconnected in this
graph. C11, wasm32 safe.

## Acceptance

`tests/test_geodesic_dijkstra.cpp`: on a unit sphere the maximum distance from
a vertex is 3.1326 (≈ π), all distances are finite, symmetric, never below
Euclidean distance, and multi-source distances dominate single-source ones.
WASM smoke checks the max distance lies in (2,4).

## Lean target

Not yet formalized.
