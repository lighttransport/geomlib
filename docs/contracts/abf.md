# Contract: angle-based flattening (ABF-style)

**Status:** implemented (Wave 3).

## Source (paper only)

- Sheffer and de Sturler, "Parameterization of Surface Meshes with Reduced
  Distortion", ACM TOMS 2001 (angle-based flattening).
- Sheffer, Levy, de Sturler, "Smooth Matching of Two Point Clouds", 2005 (ABF++
  Newton acceleration).

## Approach

Initialize with free-boundary LSCM, compute target corner angles from the 3D
mesh, then minimize the squared angle error plus a triangle-flip penalty using
per-vertex gradient descent with backtracking. Pinned vertices stay fixed.

## API / constraints

`geomlib_parameterize_abf(mesh, pin0, pin1, params, out)` and
`geomlib_abf_angle_error`. This is the nonlinear ABF objective optimized by
gradient descent, not the full ABF++ sparse Newton solver; convergence is
gradual and the function always allocates `out` itself. C11, wasm32 safe.

## Acceptance

`tests/test_abf.cpp`: a 5×5 saddle patch reduces maximum corner-angle error
from 0.034796 to 0.026488 radians in 40 iterations, with pinned UVs unchanged
and all UVs finite. WASM smoke checks a reduction.

## Lean target

Not yet formalized.
