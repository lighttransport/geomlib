# Contract: ARAP surface deformation

**Status:** implemented (Wave 3).

## Source (paper only)

- Sorkine & Alexa, "As-Rigid-As-Possible Surface Modeling", Eurographics
  Symposium on Geometry Processing 2007.

## Approach

Local step: fit a rigid rotation per vertex from rest edges to current edges
with Horn's quaternion method and the existing symmetric eigensolver. Global
step: solve an overdetermined sparse least-squares system for new positions
with CGLS, pinning handle vertices with high-weight rows. Iterate.

## API / constraints

`geomlib_deform_arap(input, handles, count, iterations, output)`; topology is
preserved, output holds deformed positions. C11, wasm32 safe.

## Acceptance

`tests/test_deform.cpp`: a 5x5 grid with one corner moved to (1.4,1.3,0.2)
reaches the target within 0.05, all vertices finite, topology unchanged.
WASM smoke deforms a plane with one handle.

## Lean target

ARAP local-step optimality (orthogonal Procrustes) is future work.
