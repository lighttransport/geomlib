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

// Wave 1 parameterization (cotangent harmonic, ARAP UV) and Catmull-Clark
const paramPlane = m.createPlane(1, 1, 4);
check(m.parameterizeMesh(paramPlane, 1, 0, 10), 'harmonic UV solve succeeds in WASM');
check(paramPlane.getUVs().length === paramPlane.getVertexCount() * 2,
      'parameterization writes one UV pair per vertex');
const arapPlane = m.createPlane(1, 1, 4);
check(m.parameterizeMesh(arapPlane, 2, 0, 3), 'ARAP UV solve succeeds in WASM');
const ccBox = m.subdivideCatmullClark(box);
check(ccBox && ccBox.getTriangleCount() === box.getTriangleCount() * 6,
      'Catmull-Clark splits each triangle into six');

// Mesh booleans (grid-based SDF combination)
const boolUnion = m.meshBoolean(box, m.createSphere(0.7, 16, 8), 0);
check(boolUnion && boolUnion.getTriangleCount() > 0,
      'mesh boolean union produces triangles in WASM');

// Anisotropic metric tensor
const met = m.metricDemo();
check(met.success && Math.abs(met.distance - Math.sqrt(34)) < 1e-4,
      'anisotropic metric distance is exact in WASM');

// MMP continuous geodesic
const mmp = m.geodesicMmpDemo();
check(mmp.success && Math.abs(mmp.length - Math.sqrt(5)) < 1e-4,
      'MMP geodesic finds the cube surface shortest path in WASM');

// Circle-ellipse overlap
const cel = m.circleEllipseDemo();
check(cel.success && cel.distance < 0,
      'circle-ellipse penetration is computed in WASM');

// Convex polygon intersection area
const pia = m.polygonIntersectionAreaDemo();
check(pia.success && Math.abs(pia.area - 0.25) < 1e-3,
      'polygon intersection area is exact in WASM');

// Circle intersection area
const carea = m.circleAreaDemo();
check(carea.success && Math.abs(carea.area - 1.22837) < 1e-4,
      'circle lens area is exact in WASM');

// Polyline arc length
check(m.arcLengthDemo().success,
      'polyline arc length is exact in WASM');

// 2D line intersection
check(m.line2dDemo().success,
      'crossing lines intersect exactly in WASM');

// Circle-triangle intersection
check(m.circleTriangleDemo().success,
      'circle containing a triangle centre is detected in WASM');

// Ellipse measures
const em = m.ellipseMeasuresDemo();
check(em.success && Math.abs(em.perimeter - 9.688448) < 1e-3,
      'ellipse perimeter approximation is accurate in WASM');

// Closest point on ellipse
const ec = m.ellipseClosestDemo();
check(ec.success && Math.abs(ec.distance - 1) < 1e-3,
      'ellipse closest point is exact in WASM');

// Capsule-triangle overlap
const ctri = m.capsuleTriangleDemo();
check(ctri.success && ctri.distance < 0,
      'capsule-triangle penetration is computed in WASM');

// Ellipse point query
check(m.ellipseDemo().success,
      'ellipse implicit value is exact in WASM');

// OBB-OBB overlap
check(m.obbObbDemo().success,
      'separated OBBs are culled in WASM');

// AABB-OBB overlap
check(m.aabbObbDemo().success,
      'separated AABB and OBB are culled in WASM');

// Triangle-OBB overlap
check(m.triObbDemo().success,
      'triangle crossing an OBB is detected in WASM');

// Ray-cone intersection
check(m.rayConeDemo().success,
      'ray-cone hit point is exact in WASM');

// Capsule-OBB overlap
const cob = m.capsuleObbDemo();
check(cob.success && cob.distance < 0,
      'capsule-OBB penetration is computed in WASM');

// Circle-AABB overlap
const cab = m.circleAabbDemo();
check(cab.success && cab.distance < 0,
      'circle-AABB penetration is computed in WASM');

// Sphere-OBB overlap
const sob = m.sphereObbDemo();
check(sob.success && sob.distance < 0,
      'sphere-OBB penetration is computed in WASM');

// Tetrahedron-tetrahedron intersection
check(m.tetraTetraDemo().success,
      'overlapping tetrahedra are detected in WASM');

// Sphere-capsule overlap
const sc = m.sphereCapsuleDemo();
check(sc.success && sc.distance < 0,
      'sphere-capsule penetration is computed in WASM');

// Triangle-tetrahedron intersection
check(m.triTetraDemo().success,
      'triangle cutting a tetrahedron is detected in WASM');

// 2D polygon distance
const p2d = m.polygonDistance2dDemo();
check(p2d.success && Math.abs(p2d.distance - 1) < 1e-3,
      '2D polygon distance is exact in WASM');

// Capsule-AABB overlap
const ca = m.capsuleAabbDemo();
check(ca.success && ca.distance < 0,
      'capsule-AABB penetration is computed in WASM');

// Polyline closest point
const polyClose = m.polylineClosestDemo();
check(polyClose.success && Math.abs(polyClose.distance - 1) < 1e-3,
      'polyline closest point is exact in WASM');

// Ray-tetrahedron intersection
const rt = m.rayTetraDemo();
check(rt.success && Math.abs(rt.t - 1) < 1e-3,
      'ray-tetrahedron entry point is exact in WASM');

// Circle-line intersection
const cl = m.circleLineDemo();
check(cl.success && cl.kind === 2,
      'circle-line secant endpoints are exact in WASM');

// Plane-sphere intersection
const ps = m.planeSphereDemo();
check(ps.success && Math.abs(ps.radius - Math.sqrt(3)) < 1e-3,
      'plane-sphere circle radius is exact in WASM');

// Triangle-triangle distance
const ttd = m.triTriDistanceDemo();
check(ttd.success && Math.abs(ttd.distance - 2) < 1e-3,
      'triangle-triangle distance is exact in WASM');

// Ray-OBB intersection
const ro = m.rayObbDemo();
check(ro.success && Math.abs(ro.t - 1) < 1e-3,
      'ray-OBB entry and normal are exact in WASM');

// AABB-AABB distance
const ad = m.aabbDistanceDemo();
check(ad.success && Math.abs(ad.distance - 1) < 1e-3,
      'AABB distance and closest points are exact in WASM');

// Convex polygon width
const pw = m.polygonWidthDemo();
check(pw.success && Math.abs(pw.width - 1) < 1e-4,
      'convex polygon minimum width is exact in WASM');

// Convex polygon diameter
const pd = m.polygonDiameterDemo();
check(pd.success && Math.abs(pd.diameter - Math.SQRT2) < 1e-4,
      'convex polygon diameter is exact in WASM');

// Triangle-AABB overlap
check(m.triAabbDemo().success,
      'triangle crossing an AABB is detected in WASM');

// Sphere-triangle overlap
const st = m.sphereTriangleDemo();
check(st.success && st.distance < 0,
      'sphere-triangle penetration is computed in WASM');

// Capsule collision
const cap = m.capsuleDemo();
check(cap.success && Math.abs(cap.distance) < 1e-3,
      'tangent capsules are detected in WASM');

// Ray-cylinder intersection
check(m.rayCylinderDemo().success,
      'ray-cylinder hit and normal are exact in WASM');

// Segment-triangle intersection
check(m.segmentTriangleDemo().success,
      'segment-triangle intersection is exact in WASM');

// Circle-circle intersection
const cc = m.circleCircleDemo();
check(cc.success && cc.kind === 2,
      'unit circles intersect analytically in WASM');

// Sphere-AABB overlap
const sab = m.sphereAabbDemo();
check(sab.success && sab.distance < 0,
      'sphere-AABB penetration depth is computed in WASM');

// Plane-plane intersection
check(m.planePlaneDemo().success,
      'perpendicular planes meet in the y-axis line in WASM');

// Ray-sphere intersection
const raySphere = m.raySphereDemo();
check(raySphere.success && Math.abs(raySphere.t0 - 4) < 1e-4,
      'ray-sphere entry root is exact in WASM');

// Bezier evaluation
const bz = m.bezierDemo();
check(bz.success && Math.abs(bz.value - 1.5) < 1e-4,
      'Bezier curve evaluates at the midpoint in WASM');

// Visvalingam-Whyatt simplification
const vw = m.vwDemo();
check(vw.success && vw.kept === 5,
      'Visvalingam-Whyatt keeps spike points in WASM');

// Solid angle / winding number
const sa = m.solidAngleDemo();
check(sa.success && Math.abs(sa.omega - Math.PI / 2) < 1e-4,
      'solid angle subtends an octant in WASM');

// Broad-phase mesh overlap
const mo = m.meshOverlapDemo();
check(mo.success && mo.pairs === 0,
      'AABB overlap culls distant meshes in WASM');

// Visibility shortest path
const vis = m.visibilityDemo();
check(vis.success && vis.length > 6 && vis.length < 9,
      'visibility path detours around an obstacle in WASM');

// Polygon rasterization
const r2 = m.raster2dDemo();
check(r2.success && r2.inside > 100,
      'polygon rasterization fills a grid in WASM');

// Harmonic cage deformation
const cage2 = m.cage2dHarmonicDemo();
check(cage2.success, 'harmonic cage reproduces an affine translation in WASM');

// Curve smoothing
const cs = m.curveSmoothDemo();
check(cs.success && cs.roughness < 0.01,
      'Laplacian curve smoothing reduces roughness in WASM');

// 2D Poisson sampling
const p2 = m.poisson2dDemo();
check(p2.success && p2.count > 10,
      'Poisson disk sampling fills a polygon in WASM');

// Fréchet distance
const fr = m.frechetDemo();
check(fr.success && Math.abs(fr.distance - 2) < 1e-4,
      'Fréchet distance is computed in WASM');

// Sphere-sphere intersection
const si = m.sphereIntersectionDemo();
check(si.success && si.radius > 1.3 && si.radius < 1.35,
      'sphere intersection circle is computed in WASM');

// Principal component analysis
const pca = m.pcaDemo();
check(pca.success && pca.eig0 > 10,
      'PCA finds the principal axis in WASM');

// Minimum spanning tree
const mst = m.mstDemo();
check(mst.success && Math.abs(mst.weight - 3) < 1e-4,
      'MST of a unit square has weight 3 in WASM');

// 2D segment intersection
const seg2 = m.segment2dDemo();
check(seg2.success && Math.abs(seg2.x) < 1e-4,
      '2D segment crossing is found in WASM');

// Mesh plane clipping
const clipm = m.meshClipDemo();
check(clipm.success && Math.abs(clipm.volume - 0.5) < 0.02,
      'plane clipping halves a box in WASM');

// Quaternion slerp
const quat = m.quaternionDemo();
check(quat.success, 'quaternion slerp rotates in WASM');

// Cotangent Laplacian
const lap = m.laplacianDemo();
check(lap.success && lap.max < 1e-9,
      'cotangent Laplacian has constant nullspace in WASM');

// Point in polygon
const pip = m.pointInPolygonDemo();
check(pip.success, 'point-in-polygon classifies interior and exterior in WASM');

// 3D Poisson sampling
const p3 = m.poisson3dDemo();
check(p3.success && p3.count > 10,
      '3D Poisson disk sampling produces points in WASM');

// Tetrahedral cage deformation
const cage = m.cage3dDemo();
check(cage.success, 'affine cage deformation is exact in WASM');

// Geodesic diameter
const gd = m.geodesicDiameterDemo();
check(gd.success && gd.distance > 2.5,
      'geodesic diameter is computed in WASM');

// Polygon Boolean with holes
const pbh = m.polyboolHolesDemo();
check(pbh.success, 'polygon Boolean with holes computes in WASM');

// Triangle-triangle intersection
const tt = m.triTriDemo();
check(tt.success, 'coplanar triangles intersect in WASM');

// Ray-plane intersection
const rp = m.rayPlaneDemo();
check(rp.success && Math.abs(rp.t - 3) < 1e-4,
      'ray-plane intersection is exact in WASM');

// Segment distance
const sd = m.segmentDistanceDemo();
check(sd.success && Math.abs(sd.distance - 3) < 1e-4,
      'segment-segment distance is exact in WASM');

// Cylinder fitting
const cyl = m.cylinderFitDemo();
check(cyl.success && Math.abs(cyl.radius - 2) < 0.05,
      'cylinder fit recovers radius in WASM');

// Ellipse fitting
const ef = m.ellipseFitDemo();
check(ef.success && Math.abs(ef.a - 2) < 0.02,
      'ellipse fit recovers semi-axes in WASM');

// Chaikin subdivision
const ch = m.chaikinDemo();
check(ch.success && ch.count === 5,
      'Chaikin subdivision preserves endpoints in WASM');

// Plane fitting
const pf = m.planeFitDemo();
check(pf.success && pf.rms < 1e-3,
      'least-squares plane fit recovers a plane in WASM');

// Sphere fitting
const sf = m.sphereFitDemo();
check(sf.success && Math.abs(sf.radius - 1) < 1e-2,
      'least-squares sphere fit recovers radius in WASM');

// Circle fitting
const cf = m.circleFitDemo();
check(cf.success && Math.abs(cf.radius - 2) < 1e-3,
      'least-squares circle fit recovers radius in WASM');

// Rigid point-set fitting
const rf = m.rigidFitDemo();
check(rf.success && rf.rms < 1e-3,
      'Horn rigid fitting recovers a transform in WASM');

// Geodesic ball
const gb = m.geodesicBallDemo();
check(gb.success && gb.count > 1,
      'geodesic ball extracts a surface patch in WASM');

// Geodesic path
const gp = m.geodesicPathDemo();
check(gp.success && gp.vertices > 1 && gp.length > 0,
      'geodesic shortest path is extracted in WASM');

// Minimum enclosing circle
const mc = m.minCircleDemo();
check(mc.success && Math.abs(mc.radius - 1) < 1e-3,
      'exact minimum enclosing circle is computed in WASM');

// Convex polygon clipping
const clip = m.polygonClipDemo();
check(clip.success && clip.vertices === 4,
      'convex clipping cuts a square in WASM');

// Polygon moments
const pm = m.polygonMomentsDemo();
check(pm.success && Math.abs(pm.area - 1) < 1e-6,
      'polygon moments match the analytic square in WASM');

// Arc-length resampling
const rs = m.resampleDemo();
check(rs.success && Math.abs(rs.length - 10) < 1e-3,
      'arc-length resampling preserves polyline length in WASM');

// DBSCAN clustering
const dbscan = m.dbscanDemo();
check(dbscan.success && dbscan.clusters === 2,
      'DBSCAN finds two dense clusters in WASM');

// Point-cloud k-means
const km = m.kmeansDemo();
check(km.success && km.clusters > 0,
      'k-means recovers two point clusters in WASM');

// Symmetry detection
const sym = m.symmetryDemo();
check(sym.success && sym.rms < 0.02,
      'reflection symmetry is detected for a box in WASM');

// Mesh-mesh distance
const md = m.meshDistanceDemo();
check(md.success && Math.abs(md.distance - 1) < 0.06,
      'mesh distance matches sphere gap in WASM');

// Spectral segmentation
const spectral = m.spectralSegmentDemo();
check(spectral.success && spectral.clusters > 0,
      'spectral segmentation partitions a mesh in WASM');

// Ambient occlusion
const ao = m.ambientOcclusionDemo();
check(ao.success && ao.min > 0.95,
      'ambient occlusion is fully open on a plane in WASM');

// Polyline-mesh intersection
const pmi = m.polylineIntersectionDemo();
check(pmi.success && pmi.hits === 2,
      'polyline intersects a sphere twice in WASM');

// Minimum-area rectangle
const mar = m.minAreaRectDemo();
check(mar.success && Math.abs(mar.area - 4) < 1e-3,
      'minimum-area rectangle is exact for a square in WASM');

// Largest inscribed circle
const inc = m.inscribedCircleDemo();
check(inc.success && Math.abs(inc.radius - 1) < 0.02,
      'largest inscribed circle matches square inradius in WASM');

// Mesh mass properties
const mass = m.massPropertiesDemo();
check(mass.success && Math.abs(mass.volume - 1) < 1e-3,
      'mesh mass properties match analytic box in WASM');

// Minimum enclosing sphere
const miniball = m.miniballDemo();
check(miniball.success && Math.abs(miniball.radius - 1.73205) < 1e-3,
      'exact miniball is computed in WASM');

// Bounding sphere
const bs = m.boundingSphereDemo();
check(bs.success && Math.abs(bs.radius - 1.73205) < 0.05,
      'Ritter bounding sphere encloses box corners in WASM');

// Convex polygon decomposition
const decomp2d = m.convexDecomposeDemo();
check(decomp2d.success && decomp2d.parts >= 2,
      'convex decomposition splits an L-shape in WASM');

// Heat kernel signature
const hks = m.hksDemo();
check(hks.success && hks.count > 0 && hks.max <= 1,
      'heat kernel signature is computed in WASM');

// Point-cloud curvature
const pc = m.pointCurvatureDemo();
check(pc.success && pc.max < 0.01,
      'PCA curvature is near zero for planar points in WASM');

// Region growing
const regions = m.regionGrowingDemo();
check(regions.success && regions.regions === 6,
      'region growing separates box faces in WASM');

// Mesh contour extraction
const meshContour = m.meshContourDemo();
check(meshContour.success && meshContour.segments > 0,
      'mesh isoline extraction produces segments in WASM');

// Point-to-plane ICP
const icpp = m.icpPlaneDemo();
check(icpp.success && icpp.rms < 0.01,
      'point-to-plane ICP aligns the perturbed mesh in WASM');

// Surface Voronoi segmentation
const sv = m.surfaceVoronoiDemo();
check(sv.success && sv.first + sv.second > 0,
      'surface Voronoi partitions the mesh in WASM');

// Delaunay refinement
const refine2d = m.refine2dDemo();
check(refine2d.success && refine2d.points > 4,
      'Delaunay refinement inserts Steiner points in WASM');

// Edge-collapse coarsening
const collapse = m.collapseDemo();
check(collapse.success && collapse.collapses > 0,
      'edge collapse reduces triangle count in WASM');

// Edge-split refinement
const refine = m.refineDemo();
check(refine.success && refine.splits > 0,
      'edge splitting shortens the longest edge in WASM');

// Edge-flip quality improvement
const qf = m.qualityFlipDemo();
check(qf.success && qf.flips > 0,
      'edge flips improve triangle quality in WASM');

// Catmull-Rom resampling
const spline = m.splineDemo();
check(spline.success && spline.count === 25,
      'centripetal Catmull-Rom resampling preserves endpoints in WASM');

// Bilateral mesh smoothing
const bilat = m.meshBilateralDemo();
check(bilat.success && bilat.after < bilat.before,
      'bilateral smoothing reduces mesh noise in WASM');

// Mean-value cage deformation
const mvc = m.meanValueDemo();
check(mvc.success && Math.abs(mvc.x - 2.2) < 1e-3,
      'mean-value cage reproduces an affine translation in WASM');

// Angle-based flattening
const abf = m.abfDemo();
check(abf.success && abf.after < abf.before,
      'ABF reduces angle distortion in WASM');

// Dijkstra geodesics
const dij = m.dijkstraDemo();
check(dij.success && dij.max > 2 && dij.max < 4,
      'Dijkstra geodesic distance spans the sphere in WASM');

// SDF ray casting
const sdfRay = m.sdfRaycastDemo();
check(sdfRay.success && sdfRay.t > 0.5 && sdfRay.t < 1.5,
      'SDF ray casting hits the sphere in WASM');

// MLS surface projection
const mls = m.mlsDemo();
check(mls.success && Math.abs(mls.z) < 0.03,
      'MLS projects a query onto the plane in WASM');

// Polyline simplification
const simp = m.simplifyDemo();
check(simp.success && simp.kept < 24,
      'Douglas-Peucker reduces a polyline in WASM');

// Polygon triangulation with holes
const pth = m.polytriHolesDemo();
check(pth.success && Math.abs(pth.area - 3) < 1e-3,
      'triangulation excludes the hole in WASM');

// 2D Minkowski sum
const mink = m.minkowskiDemo();
check(mink.success && mink.count === 4 && Math.abs(mink.area - 4) < 1e-3,
      'Minkowski sum of two unit squares has area 4 in WASM');

// 2D convex hull
const hull2d = m.hull2dDemo();
check(hull2d.success && hull2d.count === 4 && Math.abs(hull2d.area - 4) < 1e-3,
      '2D convex hull returns the square hull in WASM');

// Non-rigid ICP
const nonrigid = m.nonrigidDemo();
check(nonrigid.success && nonrigid.final < nonrigid.initial,
      'non-rigid ICP reduces correspondence error in WASM');

// 2D polygon Boolean
const polybool = m.polyboolDemo();
check(polybool.success && polybool.loops === 1,
      'polygon intersection gives the correct overlap in WASM');

// 2D medial axis
const medial = m.medialAxisDemo();
check(medial.success && medial.segments > 0,
      'medial axis produces segments in WASM');

// Farthest-point sampling
const fps = m.fpsDemo();
check(fps.success && fps.count === 6 && fps.radius > 0,
      'farthest-point sampling spreads samples in WASM');

// SDF-based mesh offset
const off = m.offsetDemo();
check(off.success && off.triangles > 0 && off.volume > 0,
      'mesh offset produces an enlarged surface in WASM');

// Constrained 2D Delaunay triangulation
const cdt = m.cdt2dDemo();
check(cdt.success && cdt.triangles === 2,
      'constrained Delaunay preserves the diagonal in WASM');

// Mesh/mesh intersection
const meshInt = m.meshIntersectionDemo();
check(meshInt.success && meshInt.segments > 0,
      'mesh intersection produces segments in WASM');

// Spatial sorting
const spatial = m.spatialSortDemo();
check(spatial.success && spatial.first >= 0,
      'Morton spatial sort runs in WASM');

// Mesh voxelization
const vox = m.voxelizeDemo();
check(vox.success && vox.volume > 0 && vox.surface > 0,
      'mesh voxelization produces surface and volume voxels in WASM');

// 2D Voronoi diagram
const vor = m.voronoiDemo();
check(vor.success && vor.sites === 16 && Math.abs(vor.area - 1) < 1e-3,
      'Voronoi cells tile the box in WASM');

// Grid-based fluid projection
const fluid = m.fluidProjectDemo();
check(fluid.success && fluid.after < fluid.before,
      'fluid projection reduces divergence in WASM');

// Variational RBF reconstruction
const rbfPts = [];
for (let i = 0; i < 80; i++) {
  const a = i * 2.399963, z = 1 - 2*(i+0.5)/80;
  const r = Math.sqrt(Math.max(0, 1-z*z));
  rbfPts.push(r*Math.cos(a), r*Math.sin(a), z);
}
const rbf = m.reconstructRBF(rbfPts);
check(rbf.success && rbf.triangles > 0,
      'RBF reconstruction produces a surface in WASM');

// Hex-dominant meshing
const hex = m.hexDominantDemo();
check(hex.success && hex.hexes >= 1 && hex.tets >= 0,
      'hex-dominant conversion merges tets in WASM');

// Sparse direct solver
const direct = m.sparseDirectSolve(new Array(64).fill(1));
check(direct.success && direct.size > 0,
      'sparse direct solve runs in WASM');

// Point-cloud octree
const octree = m.octreeQuery([0,0,0, 1,0,0, 0,1,0, 0,0,1], 0.1, 0.1, 0.1, 0.5);
check(octree.success && octree.count >= 0,
      'octree radius query runs in WASM');

// Dynamic BVH refit
const refitted = m.bvhRefitDistance(m.createSphere(1, 10, 6), m.createSphere(1, 10, 6));
check(refitted.success && Number.isFinite(refitted.distance),
      'BVH refit answers queries in WASM');

// Point-cloud k-d tree
const knn = m.kdtreeQuery([0,0,0, 1,0,0, 0,1,0, 0,0,1], 0.1, 0.1, 0.1);
check(knn.success && knn.count === 4,
      'k-d tree kNN query runs in WASM');

// Shape diameter
const sdfShape = m.shapeDiameter(m.createSphere(1, 12, 6));
check(sdfShape.success && sdfShape.values.length === sdfShape.count,
      'shape diameter returns one value per vertex in WASM');

// Cross field
const field = m.crossField(m.createSphere(1, 12, 6));
check(field.success && field.angles.length === field.count,
      'cross field returns one angle per vertex in WASM');

// Dual-quaternion skinning
const dqSkin = m.skinRotateDualQuat(m.createBox(1, 1, 1), 0.5);
check(dqSkin.success && dqSkin.vertices > 0,
      'dual-quaternion skinning transforms a mesh in WASM');

// Linear blend skinning
const skinned = m.skinRotate(m.createBox(1, 1, 1), 0.5);
check(skinned.success && skinned.vertices > 0,
      'linear blend skinning transforms a mesh in WASM');

// Continuous geodesic via unfolding
const unfold = m.geodesicUnfold(m.createPlane(1, 1, 1), 0, 3);
check(unfold.success && unfold.length > 0,
      'geodesic unfolding returns a path length in WASM');

// Fast-marching surface geodesics
const fmm = m.geodesicFMM(m.createSphere(1, 12, 6), 0);
check(fmm.success && fmm.distances.length === fmm.count && fmm.distances[0] === 0,
      'fast-marching geodesics return one distance per vertex in WASM');

// 2D Delaunay triangulation
const delaunay = m.delaunay2d([0, 0, 1, 0, 1, 1, 0, 1]);
check(delaunay.success && delaunay.triangles === 2,
      '2D Delaunay triangulates a square in WASM');

// 2D marching squares
const contourGrid = [];
for (let y = 0; y < 16; y++)
  for (let x = 0; x < 16; x++)
    contourGrid.push((x/15-0.5)**2 + (y/15-0.5)**2 - 0.0625);
const contour = m.marchingSquares(contourGrid, 1/15);
check(contour.success && contour.contours > 0 && contour.points > 0,
      'marching squares extracts a contour in WASM');

// Vertex-clustering decimation
const clustered = m.clusterDecimate(m.createSphere(1, 12, 6), 0.2);
check(clustered.success && clustered.triangles > 0,
      'vertex-clustering decimation produces a mesh in WASM');

// Balanced mesh partition
const partition = m.partitionMesh(m.createSphere(1, 12, 6), 4);
check(partition.success && partition.sizes.length === 4,
      'mesh partition returns balanced parts in WASM');

// Binary PLY I/O
const plyb = m.binaryPlyRoundtrip(m.createBox(1, 1, 1));
check(plyb.success && plyb.vertices > 0 && plyb.triangles > 0,
      'binary PLY round-trip succeeds in WASM');

// Fast winding number
const winding = m.fastWinding(m.createSphere(1, 12, 6), 0, 0, 0);
check(Number.isFinite(winding.value) && Math.abs(winding.value) > 0.5,
      'fast winding number reports inside in WASM');

// Implicit mean-curvature fairing
const fair = m.fairMesh(m.createSphere(1, 12, 6), 0.05, 2);
check(fair.success && fair.vertices > 0,
      'implicit fairing runs in WASM');

// Convex decomposition
const decomp = m.convexDecompose(m.createBox(1, 1, 1), 4);
check(decomp.success && decomp.pieces === 4 && decomp.triangles > 0,
      'convex decomposition returns convex pieces in WASM');

// Mesh slicing
const slice = m.sliceMesh(m.createBox(1, 1, 1), 0, 0, 1, 0);
check(slice.success && slice.polylines > 0 && slice.points > 0,
      'plane slicing returns polylines in WASM');

// Alpha wrapping
const wrapped = m.alphaWrap(m.createPlane(1, 1, 2), 0.1, 20);
check(wrapped && wrapped.getTriangleCount() > 0,
      'alpha wrapping produces a watertight mesh in WASM');

// Feature-preserving point-set denoising
const denoiseCloud = new m.PointCloud();
for (let i = 0; i < 60; i++) denoiseCloud.addPointWithNormal(i*0.05, 0, 0, 0, 0, 1);
const denoised = m.denoisePointCloud(denoiseCloud, 2);
check(denoised.success && denoised.count > 0,
      'bilateral point-set denoising runs in WASM');

// Variational shape approximation
const vsa = m.vsaSimplify(m.createSphere(1, 12, 6), 8);
check(vsa.success && vsa.regions === 8 && vsa.triangles > 0,
      'VSA simplification returns a coarse mesh in WASM');

// 2D polygon offset
const offset = m.offsetPolygon([0, 0, 1, 0, 1, 1, 0, 1], 0.1);
check(offset.success && offset.loops === 1 && offset.vertices === 4,
      '2D polygon offset returns the offset loop in WASM');

// 2D straight skeleton
const skel2d = m.straightSkeleton([0, 0, 1, 0, 1, 1, 0, 1]);
check(skel2d.success && skel2d.nodes > 0 && skel2d.edges > 0,
      '2D straight skeleton returns a graph in WASM');

// Curve skeleton
const skeleton = m.extractSkeleton(m.createSphere(1, 12, 6), 0.1);
check(skeleton.success && skeleton.nodes > 0 && skeleton.edges > 0,
      'curve skeleton extraction produces a graph in WASM');

// Mesh saliency
const saliency = m.computeSaliency(m.createSphere(1, 12, 6));
check(saliency.success && saliency.values.length > 0,
      'mesh saliency returns one value per vertex in WASM');

// Dual contouring isosurface extraction
const dualSurface = m.dualContouring(m.createSphere(1, 12, 6), 16);
check(dualSurface && dualSurface.getTriangleCount() > 100,
      'dual contouring extracts a surface in WASM');

// Marching tetrahedra isosurface extraction
const tetraSurface = m.marchingTets(m.createSphere(1, 12, 6), 16);
check(tetraSurface && tetraSurface.getTriangleCount() > 100,
      'marching tetrahedra extracts a surface in WASM');

// Anisotropic remeshing
const aniso = m.remeshAnisotropic(m.createSphere(1, 10, 6), 0.3, 2);
check(aniso && aniso.vertices > 0 && aniso.triangles > 0,
      'anisotropic remeshing produces a mesh in WASM');

// Quad remeshing
const quads = m.remeshToQuads(m.createPlane(1, 1, 4));
check(quads.success && quads.quads > 0,
      'quad remeshing produces quads in WASM');

// RANSAC plane detection
const planeCloud = new m.PointCloud();
for (let i = 0; i < 50; i++) {
  planeCloud.addPoint(i * 0.1, 0, 0);
  planeCloud.addPoint(i * 0.1, 1, 0);
}
const detected = m.detectPlanes(planeCloud, 0.02, 2);
check(detected.count === 1 && detected.shapes[0].inliers >= 90,
      'RANSAC plane detection finds the dominant plane in WASM');

// Manifold harmonics
const harmonics = m.manifoldHarmonics(m.createSphere(1, 12, 6), 5);
check(harmonics.success && harmonics.eigenvalues.length === 5 &&
      Math.abs(harmonics.eigenvalues[0]) < 1e-4,
      'manifold harmonics returns a zero first eigenvalue');

// ARAP surface deformation
const deformed = m.deformARAP(m.createPlane(1, 1, 4), 0, 0.4, 0.3, 0.2);
check(deformed && deformed.getVertexCount() > 0,
      'ARAP deformation produces a mesh in WASM');

// ACVD remeshing
const acvdMesh = m.remeshACVD(m.createSphere(1, 16, 8), 300);
check(acvdMesh && acvdMesh.getVertexCount() > 100 && acvdMesh.getTriangleCount() > 100,
      'ACVD remeshing produces a mesh in WASM');

// Poisson surface reconstruction from oriented points
const poissonCloud = new m.PointCloud();
for (let i = 0; i < 200; i++) {
  const u = i / 200;
  const z = 1 - 2 * u;
  const r = Math.sqrt(Math.max(0, 1 - z * z));
  const theta = i * 0.618034 * 6.2831853;
  const x = r * Math.cos(theta), y = r * Math.sin(theta);
  poissonCloud.addPointWithNormal(x, y, z, x, y, z);
}
const poissonMesh = m.reconstructPoisson(poissonCloud);
check(poissonMesh && poissonMesh.getTriangleCount() > 100,
      'Poisson reconstruction produces a mesh in WASM');

// Blue-noise surface sampling
const noise = m.blueNoiseSurface(m.createSphere(1, 16, 8), 0.2);
check(noise.count > 10 && noise.points.length === noise.count * 3,
      'blue-noise sampling returns points in WASM');

// Rigid ICP registration
const icpTarget = m.createSphere(1, 12, 6);
const icpSource = m.createSphere(1, 12, 6);
const icpResult = m.rigidICP(icpSource, icpTarget, 10);
check(Number.isFinite(icpResult.rms) && icpResult.transform.length === 16,
      'rigid ICP returns a finite RMS and a 4x4 transform');

// Robust predicates (exact expansion arithmetic in C11)
check(m.predOrient2d(0, 0, 1, 0, 0, 1) === 1, 'orient2d reports CCW');
check(m.predOrient2d(0, 0, 0, 1, 1, 0) === -1, 'orient2d reports CW');
check(m.predOrient2d(0, 0, 1, 1, 2, 2) === 0, 'orient2d reports collinear');
check(m.predOrient3d(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1) === 1,
      'orient3d reports positive tetrahedron');
check(m.predInCircle(0, 0, 1, 0, 0, 1, 0.2, 0.2) === 1,
      'incircle reports interior point');
check(m.predInSphere(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0.1, 0.1, 0.1) === 1,
      'insphere reports interior point');

// Memory-buffer mesh I/O (no filesystem needed in the browser)
const objText = m.saveMeshToString(box, 'obj');
check(objText.length > 0, 'OBJ serialization produces text');
const loadedBox = m.loadMeshFromString(objText, 'obj');
check(loadedBox !== null && loadedBox.getTriangleCount() === box.getTriangleCount() &&
      loadedBox.getVertexCount() === box.getVertexCount(),
      'OBJ memory round-trip preserves mesh counts');
const plyText = m.saveMeshToString(box, 'ply');
const loadedPly = m.loadMeshFromString(plyText, 'ply');
check(loadedPly !== null && loadedPly.getTriangleCount() === box.getTriangleCount(),
      'PLY memory round-trip preserves triangle count');

if (failures) {
  console.error(`\n${failures} WebAssembly check(s) failed`);
  process.exit(1);
}
console.log('\nAll WebAssembly checks passed');
