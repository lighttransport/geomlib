# Contract: bilateral mesh smoothing

**Status:** implemented (Wave 3).

## Source (paper only)

- Tomasi and Manduchi, "Bilateral Filtering for Gray and Color Images", ICCV
  1998; adapted to meshes via normal weighting.

## Approach

Compute vertex normals, find k nearest neighbours with a k-d tree, and update
each vertex to a Gaussian average weighted by spatial distance and normal
similarity. Updates are simultaneous (Jacobi) over the requested iterations.
The spatial width defaults to the mean edge length.

## API / constraints

`geomlib_mesh_bilateral_smooth(mesh, params)` with
`geomlib_mesh_bilateral_params_default`. Modifies vertex positions in place and
recomputes normals. k-d-tree averaging can shrink open boundaries; use it for
noise removal rather than exact boundary preservation. C11, wasm32 safe.

## Acceptance

`tests/test_mesh_bilateral.cpp`: a planar grid stays at z=0, while height noise
with RMS 0.0551 is reduced to 0.0124; all positions remain finite and normals
unit length. WASM smoke checks the noise reduction.

## Lean target

Not yet formalized.
