# Contract: surface Voronoi segmentation

**Status:** implemented (Wave 3).

## Source (paper only)

- Voronoi, 1908; the geodesic (intrinsic) restriction to a triangle mesh is a
  standard surface-partitioning primitive.

## Approach

Run one Dijkstra solve per site over the triangle-edge graph and assign each
vertex to the strictly nearest site's region (ties keep the lower index).
Region sizes are accumulated and returned with the full label array.

## API / constraints

`geomlib_surface_voronoi_compute(mesh, sites, nsites, out)` and
`geomlib_surface_voronoi_free`. Cost is O(nsites * E log V); suited to modest
site counts. C11, wasm32 safe.

## Acceptance

`tests/test_surface_voronoi.cpp`: three sites on a sphere give region sizes
65, 7 and 45; all vertices assigned, sizes sum to the vertex count, each site
owns itself, and every region is connected through mesh edges. WASM smoke
checks a two-site partition.

## Lean target

Not yet formalized.
