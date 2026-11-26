# Contract: point-cloud k-means

**Status:** implemented (Wave 3).

## Source (paper only)

- Lloyd, "Least Squares Quantization in PCM", IEEE Trans. Information Theory
  1982.

## Approach

Seed centroids with farthest-point initialization, then alternate nearest-centroid
assignment and centroid recomputation until labels stabilize or the iteration
budget is exhausted.

## API / constraints

`geomlib_pointcloud_kmeans(cloud, k, max_iterations, labels, centroids)`;
`centroids` may be NULL. Deterministic for a fixed input order. C11, wasm32
safe.

## Acceptance

`tests/test_kmeans.cpp`: two blobs at x = ±3 give clusters of 100 points each
with centroids at ±3.00, and repeated runs produce identical labels. WASM smoke
checks a 100-point two-blob cloud.

## Lean target

Not yet formalized.
