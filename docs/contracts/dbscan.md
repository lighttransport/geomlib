# Contract: DBSCAN point-cloud clustering

**Status:** implemented (Wave 3).

## Source (paper only)

- Ester, Kriegel, Sander, Xu, "A Density-Based Algorithm for Discovering
  Clusters in Large Spatial Databases with Noise", KDD 1996.

## Approach

Use a k-d tree radius query for each point's neighbourhood. Expand clusters
from dense core points through border points, marking points that never reach
the minimum neighbour count as noise. When `eps <= 0`, set it to twice the mean
nearest-neighbour distance.

## API / constraints

`geomlib_pointcloud_dbscan(cloud, params, labels, cluster_count)` writes one
label per point; noise uses `GEOMLIB_DBSCAN_NOISE`. C11, wasm32 safe.

## Acceptance

`tests/test_dbscan.cpp`: two dense blobs plus 20 sparse points give 2 clusters
with exactly 20 noise labels and all blob points clustered. WASM smoke checks
2 clusters.

## Lean target

Not yet formalized.
