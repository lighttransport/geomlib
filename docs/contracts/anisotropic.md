# Contract: anisotropic curvature-adaptive remeshing

**Status:** implemented (Wave 3).

## Sources (papers only)

- Botsch & Kobbelt, "A Remeshing Approach to Multiresolution Modeling", SGP 2004.
- Peyre & Cohen, "Surface sampling and anisotropic remeshing with adapted
  centroidal Voronoi diagrams", 2006 (metric-based remeshing concept).

## Approach

Build a per-vertex 2x2 metric tensor in the tangent plane from the existing
principal curvature module: high curvature shortens the metric length along
the corresponding principal direction. Split metric-long edges (BVH-projected
new vertices), relax vertices toward the original surface, and repeat. Split
work is bounded per pass and by a hard vertex cap (8x the input) to guarantee
termination.

## API / constraints

`geomlib_remesh_anisotropic(mesh, params)`. C11, wasm32 safe. Parameters:
target metric length, curvature scale, iterations, feature preservation.

## Acceptance

`tests/test_anisotropic.cpp`: a 153-vertex sphere refines to 1224 vertices
with all vertices finite and within 0.1 of the unit sphere. WASM smoke checks
a smaller sphere remesh.

## Lean target

Metric-tensor convergence is future work.
