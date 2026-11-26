# Contract: reflection symmetry detection

**Status:** implemented (Wave 3).

## Source (paper only)

- Principal-axis symmetry heuristics from shape analysis; reflection error via
  the closest-point query on a BVH.

## Approach

Compute the covariance of the vertices, extract its eigenvectors, and test the
three planes through the centroid whose normals are those eigenvectors. Reflect
every vertex across each plane and measure the RMS distance from the reflected
points back to the original surface; the smallest error identifies the best
symmetry plane.

## API / constraints

`geomlib_detect_symmetry(mesh, out)` returns the plane point, unit normal and
RMS error. Only principal planes are tested, so arbitrary oblique symmetries
may be missed. C11, wasm32 safe.

## Acceptance

`tests/test_symmetry.cpp`: a box gives RMS 0.00000 with a principal normal;
moving one vertex by 1.5 raises the error to 0.30844. WASM smoke checks the box.

## Lean target

Not yet formalized.
