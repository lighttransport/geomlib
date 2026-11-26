# Contract: feature-preserving point-set processing

**Status:** implemented (Wave 3).

## Sources (papers only)

- Paris, Kopf, Lang, "A Direct Approach to Bilateral Filtering of Point Sets",
  Eurographics 2007.
- Schrapp et al., "Automatic Measurement of Individual Trees" / standard
  statistical outlier removal in point-cloud processing.

## Approach

Bilateral smoothing averages neighbouring points with weights combining
spatial Gaussian distance and normal similarity, optionally estimating normals
first. Outlier removal ranks points by average k-NN distance and drops the
top `fraction`, compacting positions, normals and confidences.

## API / constraints

`geomlib_pointcloud_bilateral_smooth(cloud, params)` and
`geomlib_pointcloud_remove_outliers(cloud, fraction, k)`.
C11, wasm32 safe.

## Acceptance

`tests/test_pointset.cpp`: noisy plane points reduce mean |z| from 0.012 to
0.0006; a distant outlier is removed. WASM smoke denoises a point cloud.

## Lean target

Not yet formalized.
