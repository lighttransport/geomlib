# Contract: 2D polygon medial axis

**Status:** implemented (Wave 3).

## Source (paper only)

- Blum, "A Transformation for Extracting New Descriptors of Shape", 1967.
- The medial axis as the Voronoi diagram of a sampled polygon boundary is a
  standard computational-geometry construction.

## Approach

Sample polygon edges uniformly, build the clipped Voronoi diagram of the
samples, and treat Voronoi vertices as candidate maximal-circle centres. A
candidate is kept when it lies strictly inside the polygon and its nearest
samples belong to at least two distinct boundary edges. Kept candidates are
joined along shared Voronoi edges.

## API / constraints

`geomlib_medial_axis2d(polygon, n, spacing, out)` and `geomlib_medial2d_free`.
`spacing <= 0` selects 0.02 × bounding-box diagonal. The result is a sampling
approximation whose detail depends on spacing; it assumes a simple polygon.
C11, wasm32 safe.

## Acceptance

`tests/test_medial2d.cpp`: a unit square at spacing 0.06 gives 68 segments;
every endpoint is inside the polygon and equidistant to at least two edges,
with maximum inradius 1.000. A 2×1 rectangle gives 126 segments. WASM smoke
checks a diamond polygon.

## Lean target

Not yet formalized.
