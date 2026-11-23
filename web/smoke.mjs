// Headless smoke test for the WebAssembly module.
import createGeomlibModule from './demo/public/geomlib.js';

let failures = 0;
function check(condition, message) {
  if (condition) {
    console.log(`ok   ${message}`);
  } else {
    console.error(`FAIL ${message}`);
    failures++;
  }
}

const m = await createGeomlibModule();

// Plane + boundary loops
const plane = m.createPlane(2, 2, 4);
check(plane.getTriangleCount() > 0, 'createPlane produces triangles');
const loops = m.extractBoundaryLoops(plane);
check(loops.loopCount === 1 && loops.loops[0].vertexCount >= 4,
      'plane has one closed boundary loop');

// Closed box -> SDF -> rebuilt solid
const box = m.weldVertices(m.createBox(1, 1, 1), 1e-5);
const boxAnalysis = m.analyzeManifold(box);
check(boxAnalysis.isClosed && boxAnalysis.isManifold, 'welded box is closed and manifold');
const sdf = m.generateSDFResolution(box, 16);
check(sdf.getResolutionX() > 0 && sdf.getVoxelSize() > 0, 'SDF generation succeeds');
check(sdf.isInside(0, 0, 0), 'SDF reports the box center inside');
const rebuilt = sdf.extractSurface(0);
const rebuiltAnalysis = m.analyzeManifold(rebuilt);
check(rebuilt.getTriangleCount() > 0, 'SDF surface extraction produces triangles');
check(rebuiltAnalysis.isClosed && rebuiltAnalysis.isManifold,
      'rebuilt solid is closed and manifold');

// Geodesics (heat method)
const geodesic = m.computeGeodesicDistances(plane, 0);
check(geodesic.success && geodesic.distances.length === plane.getVertexCount(),
      'heat geodesics return one distance per vertex');
check(geodesic.distances[0] === 0, 'geodesic source distance is zero');

// Smoothing and subdivision
check(m.smoothTaubin(box, 0.5, -0.53, 2).success, 'Taubin smoothing runs');
const decimated = m.createSphere(1, 16, 8);
const beforeDecimation = decimated.getTriangleCount();
const decimation = m.decimateQEM(decimated, 0.5);
check(decimation.success && decimated.getTriangleCount() < beforeDecimation,
      'QEM decimation reduces the triangle count');
const butterfly = m.subdivideButterfly(plane, 1);
check(butterfly && butterfly.getTriangleCount() === plane.getTriangleCount() * 4,
      'modified Butterfly splits each triangle into four');
const loop = m.subdivideLoop(plane, 1);
check(loop && loop.getTriangleCount() === plane.getTriangleCount() * 4,
      'Loop subdivision splits each triangle into four');

// Reconstruction from a planar point cloud
const cloud = new m.PointCloud();
for (let y = 0; y < 4; y++) {
  for (let x = 0; x < 4; x++) {
    cloud.addPointWithNormal(x, y, 0, 0, 0, 1);
  }
}
const greedy = m.reconstructSurface(cloud, 0, 0);
check(greedy.getTriangleCount() > 0, 'greedy reconstruction produces triangles');
const marching = m.reconstructSurface(cloud, 2, 0.25);
check(marching.getTriangleCount() > 0, 'marching-cubes reconstruction produces triangles');

// Tetrahedral meshing
const tetra = m.tetrahedralizeBox(0, 0, 0, 1, 1, 1, 0.5);
check(tetra && tetra.getTetrahedronCount() > 0, 'box tetrahedralization produces cells');
const tetraStats = m.computeTetraMeshStats(tetra);
check(tetraStats.totalVolume > 0.9 && tetraStats.totalVolume < 1.1,
      'tetrahedral box volume is approximately one');

// Convex hull
const hull = m.computeConvexHull(m.createSphere(1, 16, 8));
check(hull && hull.getTriangleCount() > 0, 'convex hull produces triangles');
check(hull.getVolume() < 5.0, 'convex hull is not the bounding box');

if (failures) {
  console.error(`\n${failures} WebAssembly check(s) failed`);
  process.exit(1);
}
console.log('\nAll WebAssembly checks passed');
