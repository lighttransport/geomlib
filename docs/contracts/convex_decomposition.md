# Contract: approximate convex decomposition

**Status:** implemented (Wave 3).

## Source (paper only)

- CGAL `approximate_convex_decomposition` concept description; VSA clustering
  from Cohen-Sfady et al., SIGGRAPH 2004.

## Approach

Cluster mesh triangles with VSA planar proxies, gather each region's vertices,
and compute each region's convex hull with the existing incremental 3D hull.
The result is a set of convex triangle-mesh pieces.

## API / constraints

`geomlib_approx_convex_decomposition(mesh, region_count, out)` returning a
`geomlib_convex_decomposition` array of hull meshes. C11, wasm32 safe.

## Acceptance

`tests/test_convex_decomposition.cpp`: a unit box decomposed into four pieces
yields four hulls, at least two with positive volume. WASM smoke checks piece
and triangle counts.

## Lean target

Not yet formalized.
