# Contract: implicit mean-curvature fairing

**Status:** implemented (Wave 3).

## Source (paper only)

- Desbrun, Meyer, Schroder, Barr, "Implicit Fairing of Irregular Meshes using
  Diffusion and Curvature Flow", SIGGRAPH 1999.

## Approach

Each step solves the sparse SPD system (M + lambda L) x = M x_old for x, y and
z separately with CG; boundary vertices are pinned with high-weight identity
rows. Unconditionally stable compared with explicit Laplacian smoothing.

## API / constraints

`geomlib_fair_mesh(mesh, params)` modifies positions in place. C11, wasm32 safe.

## Acceptance

`tests/test_fairing.cpp`: noisy grid mean |z| drops from 0.0111 to 0.00039 in
3 steps and boundary vertices remain fixed. WASM smoke checks success.

## Lean target

Not yet formalized.
