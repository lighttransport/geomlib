# Contract: Shape Diameter Function

**Status:** implemented (Wave 3).

## Sources (papers only)

- Shapira, Shalom, Shamir, Cohen-Or, "Consistent Shape Maps via
  Segmented-Based Reparameterization", 2008 (shape diameter);
  Shlafman, Tal, Katz, "Symmetry of Metallic Car Body Parts", 2007.

## Approach

For each vertex, sample a cone of inward rays around the vertex normal and
take the median hit distance, normalized by the bounding-box diagonal. A
double-sided Möller-Trumbore test is used because the BVH ray tracer culls back
faces. Unused vertices have no meaningful thickness and receive zero.

## API / constraints

`geomlib_compute_shape_diameter(mesh, ray_count, out)`. C11, wasm32 safe.
Complexity is O(V * T * rays); intended for interactive-size meshes.

## Acceptance

`tests/test_shape_diameter.cpp`: a unit sphere yields normalized values
0.536–0.541 (uniform). WASM smoke checks one value per vertex.

## Lean target

Not yet formalized.
