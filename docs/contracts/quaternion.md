# Contract: quaternion operations

**Status:** implemented (Wave 3).

## Source (paper only)

- Shoemake, "Animating Rotation with Quaternion Curves", SIGGRAPH 1985 (slerp).

## Approach

Represent rotations as unit quaternions (w,x,y,z). Provide axis-angle
construction, multiplication, conjugation, normalization, 4×4 column-major
matrix conversion, shortest-path spherical interpolation and vector rotation.

## API / constraints

`geomlib_quat_*` functions in `geomlib_quaternion.h`. Slerp flips the sign of
the target quaternion when the dot product is negative to take the short arc.
C11, wasm32 safe.

## Acceptance

`tests/test_quaternion.cpp`: quaternion composition matches the product of
rotation matrices to 6e-8; slerp endpoints are exact; the half-slerp rotation
of (1,0,0) is (0.878, 0.479, 0) with unit length; q·q* is identity. WASM smoke
checks the half-slerp rotation.

## Lean target

Not yet formalized.
