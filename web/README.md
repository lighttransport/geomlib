# Geomlib WebAssembly Demo

This directory contains the WebAssembly build and Three.js demo for geomlib.

## Building

```bash
# Configure and build WebAssembly module
cd web
mkdir build && cd build
emcmake cmake ..
emmake make

# This will generate geomlib.js and geomlib.wasm in demo/public/
```

The module is built with `-sENVIRONMENT=web,node`, so the same artifact can be
exercised headlessly:

```bash
node smoke.mjs
```

The smoke test checks a plane's boundary loop, SDF generation and a closed
surface rebuilding a welded box, heat-method geodesics, Taubin smoothing,
Butterfly and Loop subdivision, point-cloud reconstruction, and box
tetrahedralization.

## Running the Demo

```bash
cd demo
bun install    # First time only
bun run dev    # Start development server

# Open http://localhost:3000/
```

## Features

### Mesh Generation
- **Sphere**: Generate UV sphere with configurable segments
- **Box**: Create box mesh
- **Cylinder**: Generate cylindrical mesh
- **Plane**: Create subdivided plane

### Point Cloud & Surface Reconstruction
- **Random Point Cloud**: Generate 500 random points in spherical distribution
- **Generate Gaussian Splats**: Simulate Gaussian splats and sample points
- **Load File (PLY/SPZ)**: Load point clouds from files
  - **PLY Format**: ASCII, binary, and compressed PLY files
  - **SPZ/Splat Format**: Gaussian splat files (compressed or uncompressed)
- **Reconstruct (Greedy)**: Tangent-plane one-ring fan triangulation with normal checks
- **Reconstruct (Marching)**: Oriented-point implicit field extracted with marching cubes
- **Reconstruct (Auto)**: Automatic algorithm selection based on point count

### Collision Proxy Generation
- **Convex Hull**: Compute convex hull of mesh
- **Simplified**: Create simplified collision mesh (50 triangles)
- **OBB**: Generate Oriented Bounding Box using PCA
- **AABB**: Generate Axis-Aligned Bounding Box

### Mesh Operations
- **Simplify 50%**: Reduce mesh complexity by 50% using quadric error metrics
- **Build BVH**: Construct Bounding Volume Hierarchy for ray tracing
- **Cast Random Rays**: Visualize ray tracing with random rays (green = hit, red = miss)
- **Rebuild from SDF**: Extract the zero level set from a generated distance field

### Visualization Options
- **Wireframe**: Toggle wireframe overlay on meshes
- **Color Components**: Color each connected component of the source mesh
- **Show Boundary Loops**: Draw closed mesh boundaries in gold
- **Show Proxy**: Display collision proxy mesh (orange)
- **Show BVH**: Visualize BVH structure with colored bounding boxes
  - **BVH Depth Slider**: Control how many levels of the BVH tree to display
  - Color-coded by depth: Red (root) → Orange → Yellow → Green → Cyan → Blue → Purple → Magenta
  - Leaf nodes are more opaque than interior nodes

### Statistics
- Mesh statistics include vertex, triangle, and edge counts, surface area, median edge length, and minimum angle. Closed meshes also show volume.
- BVH statistics (max depth, leaf count, average leaf size)

## BVH Visualization

The BVH visualization shows the spatial partitioning structure used for efficient ray tracing:

1. Click **Build BVH** to construct the hierarchy
2. Enable **Show BVH** checkbox
3. Adjust the **BVH Depth** slider to control visualization detail
   - Lower values show only top levels (coarser partitioning)
   - Higher values show deeper levels (finer partitioning)
4. Colors indicate tree depth from root (red) to leaves (varying colors)
5. Leaf nodes appear more opaque than interior nodes

This helps understand how the BVH spatially partitions the mesh for efficient ray-mesh intersection tests.

## File Loading

The demo supports loading point clouds and Gaussian splats from external files:

### Supported Formats

**PLY (Polygon File Format)**
- ASCII PLY (`.ply`)
- Binary PLY (little-endian, `.ply`)
- Compressed PLY (gzip compressed, `.ply`)

**Gaussian Splat Formats**
- SPZ (compressed splat format, `.spz`)
- Splat (uncompressed, `.splat`)

### Usage

1. Click **"Load File (PLY/SPZ)"** file input
2. Select a `.ply`, `.spz`, or `.splat` file
3. Point cloud will be loaded and visualized
4. Use reconstruction buttons to generate surface mesh

### File Format Details

**PLY Files** should contain vertex elements with at least:
- `property float x/y/z` - vertex positions
- Optional: `property float nx/ny/nz` - normals
- Optional: `property uchar red/green/blue` - colors

**SPZ/Splat Files** contain Gaussian splat data:
- Position (3 floats)
- Scale (3 floats)
- Rotation (4 floats, quaternion)
- Color (3 floats, SH or RGB)
- Opacity (1 float)

When loaded, Gaussian splats are automatically sampled to create a point cloud for surface reconstruction.

## Project Structure

```
web/
├── CMakeLists.txt       # Emscripten build configuration
├── binding.cc           # Embind JavaScript bindings
├── build/               # Build output directory
│   ├── geomlib.js
│   └── geomlib.wasm
└── demo/
    ├── package.json     # Bun/Vite project
    ├── vite.config.js   # Vite configuration
    ├── index.html       # Demo UI
    ├── style.css        # Styling
    ├── main.js          # Three.js integration
    ├── fileLoaders.js   # PLY/SPZ file parsers
    └── public/          # Static assets (wasm files)
```

## API Usage

```javascript
// The web server serves demo/public/ at its root.
import createGeomlibModule from '/geomlib.js';

const geomlib = await createGeomlibModule();

// Create mesh
const mesh = geomlib.createSphere(1.0, 32, 16);

// Get mesh data for Three.js
const vertices = mesh.getVertices();  // JavaScript array
const indices = mesh.getIndices();    // JavaScript array
const normals = mesh.getNormals();    // JavaScript array

// Generate collision proxy
const proxy = geomlib.generateCollisionProxy(mesh, 0, 0); // Convex hull

// Build BVH
const bvh = geomlib.buildBVH(mesh);

// Get BVH stats
const stats = bvh.getStats();
console.log('Max depth:', stats.maxDepth);

// Get BVH bounding boxes for visualization
const boxes = bvh.getBoundingBoxes(3); // Up to depth 3
boxes.forEach(box => {
  console.log(`Depth ${box.depth}: [${box.minX}, ${box.minY}, ${box.minZ}] to [${box.maxX}, ${box.maxY}, ${box.maxZ}]`);
});

// Trace ray
const hit = bvh.traceRay(0, 0, 5, 0, 0, -1);
if (hit.hit) {
  console.log('Hit at:', hit.x, hit.y, hit.z);
}
```

Generate a signed distance field from a closed mesh and extract a new mesh at
the desired iso-value. Keep the source mesh and its BVH alive while making
queries against that BVH.

```javascript
const sdf = geomlib.generateSDFResolution(mesh, 24);
if (sdf) {
  try {
    const rebuilt = sdf.extractSurface(0);
    if (rebuilt) {
      try {
        console.log(rebuilt.getTriangleCount());
      } finally {
        rebuilt.delete();
      }
    }
  } finally {
    sdf.delete();
  }
}
bvh.delete();
proxy.delete();
mesh.delete();
```

Call `.delete()` on owned Embind objects when finished. Topology and statistics
results are plain JavaScript objects and arrays.
