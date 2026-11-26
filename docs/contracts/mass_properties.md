# Contract: mesh mass properties

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard signed-tetrahedron integration formula for polyhedral volume,
  centre of mass and second moments.

## Approach

For each oriented triangle, integrate over the tetrahedron formed with the
origin: signed volume `dot(a, cross(b,c))/6`, first moment from its centroid,
and second moment `V/20 * (sum p_i p_j^T + (sum p)(sum p)^T)`. Translate the
inertia tensor from the origin to the centre of mass with the parallel-axis
theorem. Consistently inverted winding is normalized by flipping the signs of
all integrals.

## API / constraints

`geomlib_compute_mass_properties(mesh, out)` requires a closed mesh and unit
density. C11, wasm32 safe.

## Acceptance

`tests/test_mass_properties.cpp`: a unit box gives volume 1.0000, centre at the
origin, and inertia 0.1667 on the diagonal; a sphere gives volume 4.1013 and
inertia 1.6208 (analytic 4.1888 and 1.6755). WASM smoke checks the box.

## Lean target

Not yet formalized.
