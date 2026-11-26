# Geomlib References

This document contains academic papers, algorithms, and technical resources used in the development of geomlib.

## Mesh Healing

### Face Orientation and Topology

**Consistent Normal Orientation for Polygonal Meshes**
- Authors: Borodin et al.
- Year: 2004
- Conference: CGI 2004
- Description: Algorithms for detecting and fixing inconsistent face orientations in triangle meshes
- Used in: `geomlib_healing.c` - Face orientation propagation

**Robust Inside-Outside Segmentation using Generalized Winding Numbers**
- Authors: Jacobson et al.
- Year: 2013
- Conference: ACM Transactions on Graphics (SIGGRAPH)
- Description: Robust methods for determining inside/outside classification using winding numbers
- Used in: `geomlib_sdf.c` - generalized winding number; `geomlib_healing.c`
  uses ray casting rather than winding numbers.

**Surface Reconstruction from Unorganized Points**
- Authors: Hoppe et al.
- Year: 1992
- Conference: ACM SIGGRAPH
- Description: Classic paper on surface reconstruction and signed distance functions
- Used in: `geomlib_surface.c`, `geomlib_healing.c`

## Surface Reconstruction

**Poisson Surface Reconstruction**
- Authors: Kazhdan, Bolitho, Hoppe
- Year: 2006
- Conference: Eurographics Symposium on Geometry Processing
- Description: Implicit surface reconstruction using Poisson equation
- Influence: the oriented-point implicit field in `geomlib_surface.c` is
  conceptually related, but geomlib does not implement a Poisson solver.
  Screened Poisson reconstruction is a planned future algorithm.

**Greedy Projection Triangulation**
- Description: Fast triangulation method for point clouds
- Used in: `geomlib_surface.c` - Real-time surface reconstruction

## Remeshing

**A Remeshing Approach to Multiresolution Modeling**
- Authors: Mario Botsch, Leif Kobbelt
- Year: 2004
- Conference: Eurographics/ACM SIGGRAPH Symposium on Geometry Processing (SGP)
- Description: The foundational incremental remeshing algorithm using edge split, collapse, flip, and tangential smoothing
- Used in: `geomlib_remeshing.c` - split/collapse/flip and tangential smoothing;
  surface projection uses a BVH closest-point query (Wave 0). Adaptive sizing
  is still a local normal-variation heuristic.
- URL: https://www.pmp-library.org/remeshing.html

**Polygon Mesh Processing Library - Remeshing**
- Library: PMP (Polygon Mesh Processing)
- Description: Production-quality implementation of Botsch-Kobbelt remeshing
- Reference implementation for: Edge split/collapse thresholds, tangential smoothing, feature preservation
- URL: https://www.pmp-library.org/remeshing.html

**Adaptive Remeshing for Real-Time Mesh Deformation**
- Authors: Dunyach et al.
- Year: 2013
- Conference: Eurographics (Short Papers)
- Description: Fast adaptive remeshing for interactive applications
- URL: https://www.irit.fr/recherches/VORTEX/publications/rendu-geometrie/EGshort2013_Dunyach_et_al.pdf

**Isotropic Remeshing**
- Source: Stanford CS248 Cardinal3D
- Description: Educational resource on isotropic remeshing with clear algorithm descriptions
- Used for: Algorithm verification and parameter tuning
- URL: https://stanford-cs248.github.io/Cardinal3D/meshedit/global/remesh/

**Geometry Central - Remeshing**
- Library: Geometry Central
- Description: Modern C++ geometry processing library with remeshing implementations
- URL: https://geometry-central.net/surface/algorithms/remeshing/

**New Surface Remeshing Algorithms in CGAL**
- Year: 2025
- Description: Recent advances in CGAL including ACVD (Approximated Centroidal Voronoi Diagram) remeshing
- URL: https://www.cgal.org/2025/05/22/Surface_remeshing/

**Isotropic Remeshing with Fast and Exact Computation of Restricted Voronoi Diagram**
- Description: Efficient RVD computation for high-quality remeshing
- Used for: Understanding advanced remeshing techniques

## BVH (Bounding Volume Hierarchy)

**On fast Construction of SAH-based Bounding Volume Hierarchies**
- Authors: Wald
- Year: 2007
- Description: Surface Area Heuristic for efficient BVH construction
- Used in: `geomlib_bvh.c`

**Real-Time Ray Tracing**
- Description: General BVH traversal algorithms
- Used in: `geomlib_bvh.c` - Ray-triangle intersection and traversal

## Collision Detection

**Real-Time Collision Detection**
- Author: Christer Ericson
- Year: 2004
- Publisher: Morgan Kaufmann
- Description: Comprehensive collision detection algorithms including OBB, convex hull
- Used in: `geomlib_collision.c`

## Linear Algebra

**Numerical Recipes in C**
- Description: Classic numerical methods including Gaussian elimination, LU decomposition
- Used in: `geomlib_solver.c`

## Monte Carlo Methods

**Monte Carlo Methods in Computer Graphics**
- Description: Random sampling, importance sampling, and variance reduction
- Used in: `geomlib_montecarlo.c`, `geomlib_surface.c`

## Data Structures

**The Half-Edge Data Structure**
- Description: Standard mesh connectivity representation
- Used in: `geomlib_healing.c`, `geomlib_remeshing.c`
- Properties: Efficient one-ring queries, edge twin lookup, face traversal

**Efficient Edge-Based Data Structure**
- Description: O(n log n) twin half-edge finding using sorted edges
- Used in: `geomlib_healing.c:geomlib_halfedge_build()`

## Gaussian Splatting

**3D Gaussian Splatting for Real-Time Radiance Field Rendering**
- Authors: Kerbl et al.
- Year: 2023
- Conference: ACM SIGGRAPH
- Description: Novel 3D representation using Gaussian splats
- Used in: `geomlib_surface.c` - Point cloud generation from splats

## General References

**Physically Based Rendering: From Theory to Implementation**
- Authors: Matt Pharr, Wenzel Jakob, Greg Humphreys
- Description: Ray tracing, sampling, and geometric algorithms
- URL: https://www.pbr-book.org/

**Polygon Mesh Processing**
- Authors: Mario Botsch et al.
- Publisher: AK Peters/CRC Press
- Year: 2010
- Description: Comprehensive reference for mesh processing algorithms

## Online Resources

**libigl**
- Description: C++ geometry processing library with extensive algorithm implementations
- URL: https://libigl.github.io/

**Geometry Processing - Course Notes**
- Various universities (Stanford, MIT, etc.)
- Description: Educational materials on mesh processing, remeshing, and reconstruction

## Implementation Notes

### Remeshing Parameters (Botsch-Kobbelt)

From the original paper and production implementations:
- Edge split threshold: 4/3 × target_edge_length
- Edge collapse threshold: 4/5 × target_edge_length
- Target vertex valence: 6 (for interior vertices)
- Typical iterations: 10-15 for convergence
- Feature angle threshold: 30-45 degrees for sharp edge detection

### Quality Metrics

**Triangle Quality:**
- Ideal equilateral triangle: all angles = 60°
- Good quality: min angle > 30°, max angle < 120°
- Poor quality: min angle < 10° or max angle > 150°

**Edge Length Uniformity:**
- Uniform mesh: max_edge_length / min_edge_length < 2.0
- Adaptive mesh: Varies by curvature, typical range 0.5× to 2.0×

**Vertex Valence:**
- Interior vertices: ideal valence = 6
- Boundary vertices: ideal valence = 4
- Irregular vertices: any valence != ideal

## Algorithm Time Complexity

- Half-edge construction: O(n log n) where n = number of edges
- BVH construction: O(n log n) where n = number of triangles
- Remeshing iteration: O(n) where n = number of vertices
- Ray-BVH intersection: O(log n) where n = number of triangles
- Point cloud reconstruction (greedy): O(n²) where n = number of points
- Surface normal computation: O(n) where n = number of vertices
