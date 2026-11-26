# Contract: rigid ICP registration

**Status:** implemented (Wave 2).

## Sources (papers only)

- Besl & McKay, "A Method for Registration of 3-D Shapes", PAMI 1992.
- Horn, "Closed-form solution of absolute orientation using unit quaternions",
  JOSA 1987 — rotation fitting from the cross-covariance eigenproblem.

## Specification

Iterated closest point: for each source vertex find the closest target surface
point with the BVH, fit a rigid transform via Horn's symmetric 4x4 quaternion
matrix and the existing symmetric eigensolver, and iterate until the RMS
change is below tolerance. Parameters: max iterations, RMS tolerance, optional
correspondence-distance cutoff. Returns the final RMS and a column-major 4x4
transform.

## API / constraints

C11 double-precision internal math, float matrices at the API boundary. WASM:
yes.

## Acceptance

`tests/test_registration.cpp`: a sphere perturbed by a known rotation and
translation is realigned to RMS < 0.05; the returned transform reproduces the
alignment. WASM smoke checks the returned RMS and 4x4 matrix.

## Lean target

Rotation orthonormality and quaternion eigen-optimality are future work.
