# Contract: non-rigid ICP registration

**Status:** implemented (Wave 3).

## Source (paper only)

- Sumner, Heumann, Hornung, Roos, Seidel, "Embedded Deformation for Shape
  Manipulation", SIGGRAPH 2007.
- Besl and McKay, "A Method for Registration of 3-D Shapes", 1992 (ICP).

## Approach

Build a BVH of the target. Repeatedly find each source vertex's closest target
point, form a per-vertex displacement goal, and relax it with the mean
displacement of mesh neighbours (Jacobi blending). The fitted field is applied
to produce a new owned mesh.

## API / constraints

`geomlib_nonrigid_icp(source, target, max_iterations, smoothness, out)` with
`geomlib_nonrigid_result_free` and `geomlib_nonrigid_apply`. `smoothness` in
[0,1] controls neighbour averaging (0 = fully smooth, 1 = none). Assumes a
reasonable initial alignment and dense enough sampling; no normal or
landmark constraints. C11, wasm32 safe.

## Acceptance

`tests/test_nonrigid.cpp`: a sphere deformed by `y += 0.3 sin(x)`,
`z += 0.15x` is registered in 40 iterations, reducing closest-point RMS from
0.0747 to 0.0028 with mean vertex error 0.0984. WASM smoke checks error
reduction.

## Lean target

Not yet formalized.
