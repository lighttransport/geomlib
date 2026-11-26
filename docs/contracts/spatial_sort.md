# Contract: spatial sorting (Morton / Hilbert)

**Status:** implemented (Wave 3).

## Source (paper only)

- Morton, "A Computer Oriented Geodetic Data Base", 1966 (Z-order interleaving).
- Hilbert, "Ueber die stetige Abbildung einer Linie auf ein Flächenstück", 1891;
  the standard recursive quadrant-rotation formulation of the Hilbert order.

## Approach

Interleave coordinate bits with fixed spread/compact masks for 2D/3D Morton
codes. Hilbert 2D uses the quadrant rotation at each grid level. Sorting is a
stable qsort on (code, original index) after quantizing input bounds to the
code grid.

## API / constraints

`geomlib_morton{2,3}_encode/decode`, `geomlib_hilbert2_encode/decode`,
`geomlib_spatial_sort2` (Hilbert) and `geomlib_spatial_sort3` (Morton).
Hilbert is 2D only; 3D uses Morton. C11, wasm32 safe.

## Acceptance

`tests/test_spatial_sort.cpp`: Morton 2D/3D encode-decode round trips exact;
Hilbert on 64×64 is a bijection (4096/4096 unique codes with zero round-trip
failures); sorted consecutive points average 0.1883 apart vs 0.6406 for random
pairs. WASM smoke verifies Morton round trip and sorting.

## Lean target

Not yet formalized.
