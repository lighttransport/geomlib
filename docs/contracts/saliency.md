# Contract: mesh saliency and saliency segmentation

**Status:** implemented (Wave 3).

## Source (paper only)

- Lee, Beamish, White, "Mesh Saliency", ACM SIGGRAPH 2005.

## Approach

Estimate per-vertex mean curvature magnitude from the uniform Laplacian and
lumped triangle area; smooth it over several vertex-ring scales; accumulate
scale-weighted differences and normalize to [0,1]. Segmentation seeds regions
at the highest-saliency vertices and grows them with a deterministic priority
flood over vertex adjacency.

## API / constraints

`geomlib_compute_saliency(mesh, scales, out)` and
`geomlib_segment_by_saliency(mesh, saliency, region_count, labels)`.
C11, wasm32 safe.

## Acceptance

`tests/test_saliency.cpp`: a grid with a single raised vertex gives normalized
finite saliency (bump 0.776) and 3-region labels covering every vertex. WASM
smoke checks the returned saliency array.

## Lean target

Not yet formalized.
