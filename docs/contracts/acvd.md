# Contract: ACVD-style remeshing

**Status:** implemented (Wave 2).

## Source (paper only)

- Valette, Chassery, Prost, "Generic remeshing of 3D triangular meshes with
  metric-dependent discrete Voronoi diagrams", IEEE TVCG 2008 (ACVD).
- Du, Faber, Gunz, "Remeshing with Fast and Exact Computation of Restricted
  Voronoi Diagram" (related CVT approach).

## Approach

1. Dense area-weighted surface sampling.
2. Initialize N sites by sampling; Lloyd relaxation assigns dense samples to
   nearest sites, moves sites to cluster centroids, and projects them back to
   the input surface with the BVH.
3. Triangulate the relaxed sites with the existing greedy-projection surface
   reconstruction (tangent-plane fans with edge/normal checks).

This is an approximated CVT (no exact restricted Voronoi diagram or anisotropic
metric yet); theisotropic case is implemented. Anisotropic metric tensors are
future work.

## API / constraints

`geomlib_remesh_acvd(input, params, output)`, `params.target_vertices` in
[10, 4000], Lloyd iterations in [0, 20]. C11, wasm32 safe.

## Acceptance

`tests/test_acvd.cpp`: a 32x16 sphere remeshed to 500 sites yields 500 vertices
and area 12.5975 vs input 12.4657 (1% error). WASM smoke checks a 300-site
remesh.

## Lean target

CVT convergence proof is future work.
