# Contract: Wave 1 gap completion

**Status:** implemented.

## Algorithms

- **Texture atlas** — planarity and distortion chart segmentation (region
  growing with acceptance predicates), MaxRects best-short-side-fit and
  guillotine packing, chart rotation optimization, and a true coverage
  percentage (packed chart area / atlas extent).
- **Catmull–Clark subdivision** — face points, edge points and the valence
  vertex rule, generalized to triangular faces; each face splits into 3 quads
  triangulated into 6 triangles. Boundary vertices use the standard crease
  rule. Source: Catmull & Clark, "Recursively generated B-spline surfaces on
  arbitrary topological meshes", CAD 1978.
- **Robust repair** — `geomlib_soup_to_mesh` welds duplicate points, removes
  duplicated triangles (both windings), orients consistently;
  `geomlib_repair_nonmanifold_edges` detaches over-used edges conservatively
  and splits non-manifold vertex fans.
- **Vertex cache** — Tipsify (live-triangle FIFO with cache-stamp traversal)
  and a simple greedy cache-aware index reordering; both are deterministic.
  Sources: Kaczandmarz, "Tipsify" (2007); standard greedy cache heuristics.

## API / constraints

C11, no external dependencies; memory-buffer friendly; WASM: yes for all.

## Acceptance

- `tests/test_atlas.cpp`: all segmentation × packing combinations, coverage
  in (0,100], finite UVs.
- `tests/test_subdivision.cpp`: closed manifold preserved, 6 triangles per
  input face.
- `tests/test_repair.cpp`: soup reconstruction yields a closed oriented
  tetrahedron; over-used edge detached and non-manifold elements removed.
- `tests/test_optimize.cpp`: all three vcache algorithms preserve the triangle
  set.

## Lean target

Marching-cubes and predicate layers already formalized; subdivision mask
identities and atlas packing optimality are future work.
