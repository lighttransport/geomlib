# Contract: anisotropic metric tensors

**Status:** implemented (Wave 3).

## Source (paper only)

- Valette, Chassery and Prost, "Generic remeshing of 3D triangular meshes with
  metric-dependent discrete Voronoi diagrams", IEEE TVCG 2008 (metric tensors
  driving anisotropic CVT).

## Approach

Represent a metric as a symmetric positive-definite 3×3 matrix. Provide
identity/diagonal/frame construction, the induced norm and distance,
determinant and anisotropy ratio, affine interpolation in the SPD cone, and a
curvature-driven surface metric: tangent sizes shrink with principal curvature
and the normal direction receives its own scale. A per-vertex mesh field is
built from the existing curvature module.

## API / constraints

`geomlib_metric3` and the `geomlib_metric_*` helpers use double precision.
Curvature construction clamps tangent sizes to at least 10% of the target size.
C11, wasm32 safe.

## Acceptance

`tests/test_metric.cpp`: the identity metric is Euclidean (7); diag(4,1,0.25)
gives sqrt(34) with determinant 1 and anisotropy 4; frame rotation preserves
distance; a curvature metric is anisotropic; a sphere metric field has positive
determinant. WASM smoke checks the weighted distance.

## Lean target

Not yet formalized.
