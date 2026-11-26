# Contract: dual-quaternion skinning

**Status:** implemented (Wave 3).

## Source (paper only)

- Kavan, Skala, Slavik, Zara, "Skinning with Dual Quaternions", SIGGRAPH 2007.

## Approach

Convert each rigid bone transform to a unit dual quaternion (real rotation,
dual translation), blend both parts with vertex weights, normalize the blended
real part, and transform each rest vertex by the resulting screw motion.

## API / constraints

`geomlib_skin_mesh_dualquat(input, weights, bone_transforms, bone_count, output)`
using the same `geomlib_skin_weight` structure as LBS. C11, wasm32 safe.

## Acceptance

`tests/test_dualquat_skinning.cpp`: a single 90-degree bone rotates exactly;
a 50/50 identity/rotation blend preserves the point's radius (1.0 vs 0.7071
for LBS). WASM smoke rotates a box.

## Lean target

Not yet formalized.
