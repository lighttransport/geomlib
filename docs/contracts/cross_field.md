# Contract: cross fields

**Status:** implemented (Wave 3).

## Source (paper only)

- Ray & Sokolov, "Robust 4-RoSy Field Generation on Surface Meshes", 2014.

## Approach

Seed a 4-directional field from principal curvature directions projected into
each vertex's tangent frame, then relax angles across edges under the 90-degree
symmetry quotient. Frame-transition corrections between tangent planes are not
yet included, so neighbouring angles are consistent modulo 90 degrees rather
than numerically identical.

## API / constraints

`geomlib_compute_cross_field(mesh, out)` returning per-vertex angles and local
frame hints. C11, wasm32 safe.

## Acceptance

`tests/test_cross_field.cpp`: sphere field is finite and respects the
90-degree quotient. WASM smoke checks one angle per vertex.

## Lean target

Not yet formalized.
