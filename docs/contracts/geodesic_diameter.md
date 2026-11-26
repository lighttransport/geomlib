# Contract: geodesic diameter

**Status:** implemented (Wave 3).

## Source (paper only)

- Graph eccentricity on triangle meshes; double-sweep heuristic from Eppstein,
  "Square Rooting with the Squeaky Wheel" (2000) style farthest-point sweeps.

## Approach

`geomlib_geodesic_diameter_exact` runs Dijkstra from every vertex and takes the
largest distance. `geomlib_geodesic_diameter_sweep` runs two Dijkstra sweeps
(from an arbitrary vertex, then from its farthest vertex), which is much faster
but not guaranteed exact on general graphs.

## API / constraints

Both return the endpoint pair and distance. Exact search is O(V · E log V);
suited to small meshes. C11, wasm32 safe.

## Acceptance

`tests/test_geodesic_diameter.cpp`: a sphere gives exact diameter 3.1153 (near
π), the double sweep agrees exactly, endpoints are valid, and the distance
matches a fresh Dijkstra query. WASM smoke checks a value above 2.5.

## Lean target

Not yet formalized.
