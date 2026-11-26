# Contract: point-cloud local curvature

**Status:** implemented (Wave 3).

## Source (paper only)

- Hoppe, DeRose, Duchamp, McDonald, Stuetzle, "Surface Reconstruction from
  Unorganized Points", SIGGRAPH 1992 (local PCA normals).

## Approach

For each point, gather k nearest neighbours with a k-d tree, compute the
centred 3×3 covariance, and use a cyclic Jacobi eigen decomposition. The
normalized smallest eigenvalue is the curvature and its eigenvector is the
local normal (sign is arbitrary).

## API / constraints

`geomlib_pointcloud_curvature(cloud, k, curvature, normals)`; either output
array may be NULL. `k <= 0` selects 16 neighbours. C11, wasm32 safe.

## Acceptance

`tests/test_point_curvature.cpp`: a planar cloud gives max curvature 0.000000
and vertical normals; a sphere gives average curvature 0.14830 and normals
aligned with the radial direction (|dot| = 1.000). WASM smoke checks the plane.

## Lean target

Not yet formalized.
