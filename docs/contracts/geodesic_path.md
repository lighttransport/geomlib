# Contract: geodesic shortest-path polyline

**Status:** implemented (Wave 3).

## Source (paper only)

- Dijkstra, 1959; graph shortest paths on triangle meshes.

## Approach

Build the undirected triangle-edge graph with Euclidean edge lengths, run
Dijkstra with predecessor tracking, then reconstruct the vertex chain from
target back to source. The result is the exact edge-graph shortest path, not
the continuous (MMP) geodesic.

## API / constraints

`geomlib_geodesic_path_between(mesh, source, target, out)` with
`geomlib_geodesic_path_free`. A path between a vertex and itself has one
vertex and zero length. C11, wasm32 safe.

## Acceptance

`tests/test_geodesic_path.cpp`: a sphere path from vertex 0 to the midpoint
vertex has 11 vertices, length 1.5607, endpoints correct, consecutive vertices
connected by mesh edges, and length matching Dijkstra. WASM smoke checks a
multi-vertex path.

## Lean target

Not yet formalized.
