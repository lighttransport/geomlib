# Contract: robust geometric predicates

**Status:** implemented (Wave 0). **Source:** Jonathan R. Shewchuk,
"Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
Predicates", Discrete & Computational Geometry 18(3):305–363, 1997. Paper only.

## Specification

- `geomlib_pred_orient2d/3d` return the exact sign of the orientation
  determinant for finite IEEE-754 double inputs.
- `geomlib_pred_incircle` / `geomlib_pred_insphere` return the exact sign of
  the in-circle / in-sphere determinant, with "inside is positive" for a
  positively oriented reference simplex.
- `geomlib_pred_det_sign` returns the exact sign of any n×n determinant,
  1 ≤ n ≤ 5.

## Derivation

A conservative floating-point filter (`32 · DBL_EPSILON · Σ|terms|`) decides
most cases. The exact fallback uses error-free transformations (`two_sum`,
Dekker `two_product`) and recursive Laplace expansion in non-overlapping
double expansions. No target library code or predicate tables were copied.

## API / constraints

C11 only; `int64_t`/`uint64_t` are available for future fixed-width extensions;
no `long double`, no compiler extensions, no deep recursion (depth ≤ 5).
WASM: yes.

## Acceptance

`tests/test_predicates.cpp` checks 480 golden signs computed independently
with Python `fractions.Fraction` (`scripts/gen_pred_vectors.py`), covering
large near-degenerate integer coordinates.

## Lean target

`formal/GeomLibFormal/Predicates.lean` specifies the integer determinants and
verifies the equivalence of the homogeneous and cross-product formulations on
a finite domain. General algebraic identities are future work.
