# Contract: spectral mesh segmentation

**Status:** implemented (Wave 3).

## Source (paper only)

- Shi and Malik, "Normalized Cuts and Image Segmentation", CVPR 1997.
- Ng, Jordan, Weiss, "On Spectral Clustering: Analysis and an Algorithm",
  NIPS 2001.

## Approach

Embed each vertex with the lowest nontrivial Laplace-Beltrami eigenvectors
(skipping the constant one), seed k-means with the largest-norm embedding
vectors, and iterate assignments and centroid updates deterministically.

## API / constraints

`geomlib_spectral_segment(mesh, k, eigen_count, labels)` writes one label per
vertex; `eigen_count` must be at least k and at least 2. Requires the harmonics
solver to converge. C11, wasm32 safe.

## Acceptance

`tests/test_spectral_segment.cpp`: a 7×7 planar grid split into 2 clusters
gives sizes 9 and 40 with mean x 1.000 and 3.450. WASM smoke checks both
clusters are non-empty.

## Lean target

Not yet formalized.
