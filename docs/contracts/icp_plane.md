# Contract: point-to-plane ICP registration

**Status:** implemented (Wave 3).

## Source (paper only)

- Chen and Medioni, "Object Modelling by Registration of Multiple Range
  Images", Image and Vision Computing 1992.

## Approach

At each iteration, find closest target points and normals for the transformed
source vertices, then solve the linearized rigid fit in the six parameters
(rotation vector, translation) by minimizing normal distances. The incremental
transform is composed with the accumulated one via Rodrigues' rotation formula.

## API / constraints

`geomlib_icp_point_to_plane(source, target, params, transform, final_rms)` with
`geomlib_icp_plane_params_default`. Returns a column-major `geomlib_mat4`.
Assumes a reasonable initial alignment and reliable target normals. C11, wasm32
safe.

## Acceptance

`tests/test_icp_plane.cpp`: a sphere translated by (0.15, -0.10, 0.05) is
realigned from RMS 0.11116 to 0.00000; the reported RMS matches direct
measurement. WASM smoke checks RMS below 0.01.

## Lean target

Not yet formalized.
