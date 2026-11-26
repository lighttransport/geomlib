# Contract: linear blend skinning

**Status:** implemented (Wave 3).

## Source (paper only)

- Magnenat-Thalmann et al., "Joint-dependent local deformations for hand
  animation and garment animation", 1988 (classic LBS).

## Approach

Each vertex stores up to four bone indices and weights; the deformed position
is the normalized weighted blend of the corresponding rigid transforms applied
to the rest position.

## API / constraints

`geomlib_skin_mesh(input, weights, bone_transforms, bone_count, output)`.
C11, wasm32 safe.

## Acceptance

`tests/test_skinning.cpp`: a single 90-degree bone rotates vertices exactly;
two blended transforms produce the midpoint. WASM smoke rotates a box.

## Lean target

Not yet formalized.
