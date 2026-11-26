# Contract: grid-based incompressible fluid solver

**Status:** implemented (Wave 3).

## Source (paper only)

- Stam, "Stable Fluids", SIGGRAPH 1999.

## Approach

Collocated velocity grid with semi-Lagrangian advection (trilinear
backtracing) and a pressure-projection step: forward-difference divergence,
a 7-point Poisson solve via RCM + banded Cholesky, and a backward-difference
pressure gradient so `D_forward · G_backward` equals the Poisson stencil.
The Neumann nullspace is regularized and the RHS mean removed.

## API / constraints

`geomlib_fluid_create/free`, `geomlib_fluid_project`, `geomlib_fluid_advect`,
`geomlib_fluid_rms_divergence`. Grid dimensions are caller-chosen; the Poisson
solve cost grows with grid size. C11, wasm32 safe.

## Acceptance

`tests/test_fluid.cpp`: on an 8³ field the interior divergence RMS drops from
6.47 to 0.21 (>90% reduction) and advection keeps all velocities finite.
WASM smoke runs the same projection. Boundary cells of the collocated grid
retain some divergence — a known limitation of this discretization; a MAC
staggered grid would be needed for exact boundary incompressibility.

## Lean target

Not yet formalized.
