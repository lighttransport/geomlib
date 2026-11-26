# Contract: parameterization completion

**Status:** implemented (Wave 1).

## Sources (papers only)

- Tutte, "How to draw a graph", Proc. LMS 1963 — harmonic/Tutte embedding.
- Levy, Petitjean, Ray, Maillot, "Least Squares Conformal Maps for Automatic
  Texture Atlas Generation", SIGGRAPH 2002 — free-boundary LSCM.
- Liu, Zhang, Gotsman, Gu, "A Local/Global Approach to Mesh Parameterization",
  SGP 2008 — ARAP UV.

## Specification

- `geomlib_parameterize_harmonic`: cotangent stiffness solve with Dirichlet
  boundary pinned to a circle/square, solved with the sparse CG solver.
- `geomlib_parameterize_lscm_free`: one complex conformality equation per
  triangle, solved as a real sparse least-squares system with two pinned
  vertices (removed translation/rotation/scale nullspace).
- `geomlib_parameterize_arap`: local best-fit rotations and a global sparse
  least-squares solve, initialized from the harmonic solution with a uniform
  local tangent frame per vertex and rest edges scaled to the [0,1] atlas.
- Distortion statistics now report angle RMS, L2/Linf stretch and area RMS.

## API / constraints

C11, double precision, uses the Wave 0 sparse module. WASM: yes.

## Acceptance

`tests/test_parameterization.cpp`: a triangulated disk, no flipped triangles,
exact boundary pins, correct pin positions, flip detection; WASM smoke covers
harmonic and ARAP.

## Lean target

Discrete conformal-energy identities are future work; L2 in-circle axioms are
formalized in `formal/GeomLibFormal/Delaunay.lean`.
