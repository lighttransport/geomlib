# Contract: geodesic ball extraction

**Status:** implemented (Wave 3).

## Source (paper only)

- Dijkstra, 1959; intrinsic geodesic disks on triangle meshes.

## Approach

Run single-source Dijkstra on the triangle-edge graph and include every vertex
whose distance to the source is at most the requested radius. Returns both a
mask and the list of included indices.

## API / constraints

`geomlib_extract_geodesic_ball(mesh, source, radius, out)` with
`geomlib_geodesic_ball_free`. Radius zero includes only the source (unless
zero-length edges exist, as at sphere poles). C11, wasm32 safe.

## Acceptance

`tests/test_geodesic_ball.cpp`: on a 5×5 grid, radii 0/0.5/1/2 include 1/1/3/6
vertices and radius 10 includes all 25. WASM smoke checks a radius-2 patch.

## Lean target

Not yet formalized.
