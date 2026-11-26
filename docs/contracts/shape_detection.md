# Contract: RANSAC shape detection

**Status:** implemented (Wave 3).

## Source (paper only)

- Fischler & Bolles, "Random Sample Consensus: A Paradigm for Model Fitting
  with Applications to Image Analysis and Automated Cartography", CACM 1981.
- Schnabel et al., "Efficient RANSAC for Point-Cloud Shape Detection", CGF 2007
  (related efficient shape detection framework).

## Approach

For each shape type (plane, sphere, cylinder), draw minimal random samples
from the remaining points, build a primitive, count inliers within a distance
threshold, keep the best consensus set over 300 iterations, refine the model
on the inliers (PCA for planes, algebraic fit for spheres), remove the inliers
and repeat up to `max_shapes`.

## API / constraints

`geomlib_detect_shapes_ransac(cloud, type, threshold, max_shapes, seed, out_count)`.
Deterministic for a fixed seed. C11, wasm32 safe.

## Acceptance

`tests/test_shape_detection.cpp`: 200 points on two parallel planes at z=0 and
z=1 yield two detected planes with 100 inliers each, vertical normals and unit
separation. WASM smoke detects a dominant plane.

## Lean target

Not applicable (statistical estimation).
