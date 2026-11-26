# Contract: visibility-graph shortest path

**Status:** implemented (Wave 3).

## Source (paper only)

- de Berg et al., "Computational Geometry: Algorithms and Applications",
  visibility graph for shortest paths amid polygonal obstacles.

## Approach

Build a graph whose nodes are the start, goal and all obstacle vertices. Connect
mutually visible pairs (segment does not properly cross any obstacle edge), then
run Dijkstra. Travelling along an obstacle edge is permitted.

## API / constraints

`geomlib_visibility_shortest_path(start, goal, obstacles, sizes, nh, out)` for
simple polygon obstacles. O(V²) visibility tests. C11, wasm32 safe.

## Acceptance

`tests/test_visibility.cpp`: an empty field gives distance 6.0000; a unit square
between (−3,0) and (3,0) gives a 4-node detour of length 6.4721. WASM smoke
checks a detour length in (6,9).

## Lean target

Not yet formalized.
