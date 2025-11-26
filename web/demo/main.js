import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import createGeomlibModule from './public/geomlib.js';
import { loadPointCloudFile } from './fileLoaders.js';

// ============================================================================
// Global State
// ============================================================================

let scene, camera, renderer, controls;
let geomlib;
let currentMesh = null;
let currentProxy = null;
let currentBVH = null;
let currentPointCloud = null;
let mainMeshGroup = null;
let proxyMeshGroup = null;
let rayLines = null;
let bvhBoxesGroup = null;
let pointCloudPoints = null;

const state = {
  showWireframe: true,
  showProxy: false,
  showBVH: false,
  showPoints: false,
  bvhMaxDepth: 3
};

// ============================================================================
// Three.js Setup
// ============================================================================

function initThreeJS() {
  // Scene
  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x1a1a1a);

  // Camera
  camera = new THREE.PerspectiveCamera(
    60,
    window.innerWidth / window.innerHeight,
    0.1,
    1000
  );
  camera.position.set(3, 3, 3);
  camera.lookAt(0, 0, 0);

  // Renderer
  const container = document.getElementById('canvas-container');
  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(window.devicePixelRatio);
  container.appendChild(renderer.domElement);

  // Controls
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.05;
  controls.minDistance = 1;
  controls.maxDistance = 20;

  // Lights
  const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
  scene.add(ambientLight);

  const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
  directionalLight.position.set(5, 10, 5);
  scene.add(directionalLight);

  const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.4);
  directionalLight2.position.set(-5, 5, -5);
  scene.add(directionalLight2);

  // Grid
  const grid = new THREE.GridHelper(10, 10, 0x444444, 0x222222);
  scene.add(grid);

  // Axes helper
  const axes = new THREE.AxesHelper(2);
  scene.add(axes);

  // Handle window resize
  window.addEventListener('resize', onWindowResize);

  // Animation loop
  animate();
}

function onWindowResize() {
  const container = document.getElementById('canvas-container');
  camera.aspect = container.clientWidth / container.clientHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(container.clientWidth, container.clientHeight);
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}

// ============================================================================
// Geomlib Mesh Visualization
// ============================================================================

function createThreeMesh(geomlibMesh, color = 0x4CAF50, wireframe = false) {
  const vertices = geomlibMesh.getVertices();
  const indices = geomlibMesh.getIndices();
  const normals = geomlibMesh.getNormals();

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
  geometry.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3));
  geometry.setIndex(Array.from(indices));

  const material = new THREE.MeshStandardMaterial({
    color: color,
    wireframe: wireframe,
    side: THREE.DoubleSide,
    flatShading: false,
    metalness: 0.3,
    roughness: 0.7
  });

  return new THREE.Mesh(geometry, material);
}

function displayMesh(geomlibMesh, isProxy = false) {
  // Remove old mesh
  if (!isProxy && mainMeshGroup) {
    scene.remove(mainMeshGroup);
    mainMeshGroup = null;
  } else if (isProxy && proxyMeshGroup) {
    scene.remove(proxyMeshGroup);
    proxyMeshGroup = null;
  }

  // Create mesh group
  const group = new THREE.Group();

  // Solid mesh
  const solidMesh = createThreeMesh(geomlibMesh, isProxy ? 0xff9800 : 0x4CAF50, false);
  group.add(solidMesh);

  // Wireframe overlay
  if (state.showWireframe) {
    const wireframeMesh = createThreeMesh(geomlibMesh, 0x000000, true);
    wireframeMesh.material.opacity = 0.3;
    wireframeMesh.material.transparent = true;
    group.add(wireframeMesh);
  }

  scene.add(group);

  if (isProxy) {
    proxyMeshGroup = group;
    proxyMeshGroup.visible = state.showProxy;
  } else {
    mainMeshGroup = group;
  }

  // Update statistics
  updateMeshStats(geomlibMesh);
}

function updateMeshStats(geomlibMesh) {
  const stats = document.getElementById('mesh-stats');
  const verts = geomlibMesh.getVertexCount();
  const tris = geomlibMesh.getTriangleCount();
  const area = geomlibMesh.getSurfaceArea();
  const volume = geomlibMesh.getVolume();

  stats.innerHTML = `
    <strong>Mesh:</strong><br>
    Vertices: ${verts}<br>
    Triangles: ${tris}<br>
    Surface Area: ${area.toFixed(3)}<br>
    Volume: ${volume.toFixed(3)}
  `;
}

function displayHealingStats(issues) {
  const statsDiv = document.getElementById('healing-stats');
  const statusColor = (isGood) => isGood ? '#4CAF50' : '#ff6b6b';

  statsDiv.innerHTML = `
    <strong>Mesh Topology:</strong><br>
    Degenerate Triangles: ${issues.degenerateTriangles}<br>
    Duplicate Vertices: ${issues.duplicateVertices}<br>
    Boundary Edges: ${issues.boundaryEdges}<br>
    Non-manifold Edges: ${issues.nonManifoldEdges}<br>
    Non-manifold Vertices: ${issues.nonManifoldVertices}<br>
    Inconsistent Faces: ${issues.inconsistentFaces}<br>
    <span style="color: ${statusColor(issues.isManifold)}">
      ${issues.isManifold ? '✓' : '✗'} Manifold
    </span><br>
    <span style="color: ${statusColor(issues.isWatertight)}">
      ${issues.isWatertight ? '✓' : '✗'} Watertight
    </span><br>
    <span style="color: ${statusColor(issues.isOriented)}">
      ${issues.isOriented ? '✓' : '✗'} Oriented
    </span>
  `;
}

function displayRemeshingStats(stats) {
  const statsDiv = document.getElementById('remeshing-stats');

  statsDiv.innerHTML = `
    <strong>Remeshing Quality:</strong><br>
    Mean Edge Length: ${stats.meanEdgeLength.toFixed(4)}<br>
    Min Edge Length: ${stats.minEdgeLength.toFixed(4)}<br>
    Max Edge Length: ${stats.maxEdgeLength.toFixed(4)}<br>
    Min Angle: ${stats.minAngle.toFixed(1)}°<br>
    Max Angle: ${stats.maxAngle.toFixed(1)}°<br>
    Irregular Vertices: ${stats.irregularVertices}
  `;
}

// ============================================================================
// UI Event Handlers
// ============================================================================

function setupEventHandlers() {
  // Mesh generation
  document.getElementById('btn-sphere').addEventListener('click', () => {
    currentMesh = geomlib.createSphere(1.0, 32, 16);
    displayMesh(currentMesh);
    clearProxy();
    clearBVH();
  });

  document.getElementById('btn-box').addEventListener('click', () => {
    currentMesh = geomlib.createBox(1.5, 1.5, 1.5);
    displayMesh(currentMesh);
    clearProxy();
    clearBVH();
  });

  document.getElementById('btn-cylinder').addEventListener('click', () => {
    currentMesh = geomlib.createCylinder(0.8, 2.0, 32);
    displayMesh(currentMesh);
    clearProxy();
    clearBVH();
  });

  document.getElementById('btn-plane').addEventListener('click', () => {
    currentMesh = geomlib.createPlane(3.0, 3.0, 10);
    displayMesh(currentMesh);
    clearProxy();
    clearBVH();
    clearPointCloud();
  });

  // Point cloud and surface reconstruction
  document.getElementById('btn-random-cloud').addEventListener('click', () => {
    // Generate random point cloud in spherical distribution
    const pc = new geomlib.PointCloud();
    const count = 500;
    for (let i = 0; i < count; i++) {
      const theta = Math.random() * Math.PI * 2;
      const phi = Math.random() * Math.PI;
      const r = 0.8 + Math.random() * 0.4;

      const x = r * Math.sin(phi) * Math.cos(theta);
      const y = r * Math.sin(phi) * Math.sin(theta);
      const z = r * Math.cos(phi);

      pc.addPoint(x, y, z);
    }

    currentPointCloud = pc;
    visualizePointCloud(currentPointCloud);
    state.showPoints = true;
    document.getElementById('chk-points').checked = true;
    updatePointCloudVisibility();
    clearMesh();
    console.log('Generated random point cloud with', count, 'points');
  });

  document.getElementById('btn-gaussian-splats').addEventListener('click', () => {
    // Generate Gaussian splat positions (simulated)
    const count = 200;
    const positions = [];

    for (let i = 0; i < count; i++) {
      const theta = Math.random() * Math.PI * 2;
      const phi = Math.random() * Math.PI;
      const r = 0.8 + Math.random() * 0.4;

      const x = r * Math.sin(phi) * Math.cos(theta);
      const y = r * Math.sin(phi) * Math.sin(theta);
      const z = r * Math.cos(phi);

      positions.push(x, y, z);
    }

    // Sample points from Gaussians
    currentPointCloud = geomlib.gaussiansSamplePoints(positions, 5);
    visualizePointCloud(currentPointCloud);
    state.showPoints = true;
    document.getElementById('chk-points').checked = true;
    updatePointCloudVisibility();
    clearMesh();
    console.log('Generated', currentPointCloud.getPointCount(), 'points from Gaussian splats');
  });

  document.getElementById('btn-reconstruct-greedy').addEventListener('click', () => {
    if (!currentPointCloud) {
      alert('Create a point cloud first!');
      return;
    }

    console.log('Reconstructing surface with greedy projection...');
    const mesh = geomlib.reconstructSurface(currentPointCloud, 0, 0.0); // GREEDY_PROJECTION
    currentMesh = mesh;
    displayMesh(currentMesh);
    console.log('Surface reconstructed:', currentMesh.getTriangleCount(), 'triangles');
  });

  document.getElementById('btn-reconstruct-marching').addEventListener('click', () => {
    if (!currentPointCloud) {
      alert('Create a point cloud first!');
      return;
    }

    console.log('Reconstructing surface with marching cubes...');
    const mesh = geomlib.reconstructSurface(currentPointCloud, 2, 0.15); // MARCHING_CUBES, voxel size
    currentMesh = mesh;
    displayMesh(currentMesh);
    console.log('Surface reconstructed:', currentMesh.getTriangleCount(), 'triangles');
  });

  document.getElementById('btn-reconstruct-auto').addEventListener('click', () => {
    if (!currentPointCloud) {
      alert('Create a point cloud first!');
      return;
    }

    console.log('Reconstructing surface with auto settings...');
    const mesh = geomlib.reconstructSurfaceAuto(currentPointCloud, 0);
    currentMesh = mesh;
    displayMesh(currentMesh);
    console.log('Surface reconstructed:', currentMesh.getTriangleCount(), 'triangles');
  });

  // Collision proxy generation
  document.getElementById('btn-convex').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }
    currentProxy = geomlib.generateCollisionProxy(currentMesh, 0, 0); // CONVEX_HULL
    displayMesh(currentProxy, true);
    state.showProxy = true;
    updateProxyVisibility();
  });

  document.getElementById('btn-simplified').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }
    currentProxy = geomlib.generateCollisionProxy(currentMesh, 1, 50); // SIMPLIFIED
    displayMesh(currentProxy, true);
    state.showProxy = true;
    updateProxyVisibility();
  });

  document.getElementById('btn-obb').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }
    currentProxy = geomlib.generateCollisionProxy(currentMesh, 2, 0); // OBB
    displayMesh(currentProxy, true);
    state.showProxy = true;
    updateProxyVisibility();
  });

  document.getElementById('btn-aabb').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }
    currentProxy = geomlib.generateCollisionProxy(currentMesh, 3, 0); // AABB
    displayMesh(currentProxy, true);
    state.showProxy = true;
    updateProxyVisibility();
  });

  // Mesh operations
  document.getElementById('btn-simplify').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }
    const simplified = geomlib.simplifyMesh(currentMesh, 0.5);
    currentMesh = simplified;
    displayMesh(currentMesh);
    clearProxy();
    clearBVH();
  });

  document.getElementById('btn-bvh').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    console.log('Building BVH...');
    currentBVH = geomlib.buildBVH(currentMesh);

    const stats = currentBVH.getStats();
    const statsDiv = document.getElementById('bvh-stats');
    statsDiv.innerHTML = `
      <strong>BVH:</strong><br>
      Max Depth: ${stats.maxDepth}<br>
      Leaf Count: ${stats.leafCount}<br>
      Avg Leaf Size: ${stats.avgLeafSize.toFixed(2)}
    `;

    // Enable BVH visualization controls
    document.getElementById('chk-bvh').disabled = false;

    // Update depth slider max to actual BVH depth
    const depthSlider = document.getElementById('bvh-depth-slider');
    depthSlider.max = stats.maxDepth;
    depthSlider.value = Math.min(3, stats.maxDepth);
    state.bvhMaxDepth = parseInt(depthSlider.value);
    document.getElementById('bvh-depth-value').textContent = state.bvhMaxDepth;

    console.log('BVH built successfully', stats);
  });

  document.getElementById('btn-ray').addEventListener('click', () => {
    if (!currentBVH) {
      alert('Build BVH first!');
      return;
    }

    castRandomRays(50);
  });

  // Mesh healing
  document.getElementById('btn-analyze').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    const issues = geomlib.analyzeMeshTopology(currentMesh);
    displayHealingStats(issues);
    console.log('Mesh topology analysis:', issues);
  });

  document.getElementById('btn-fix-orientation').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    geomlib.fixFaceOrientation(currentMesh, 0); // 0 = PROPAGATION strategy
    displayMesh(currentMesh);
    console.log('Fixed face orientation');
  });

  document.getElementById('btn-remove-degenerates').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    const removed = geomlib.removeDegenerateTriangles(currentMesh);
    displayMesh(currentMesh);
    console.log('Removed', removed, 'degenerate triangles');
  });

  document.getElementById('btn-remove-duplicates').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    const removed = geomlib.removeDuplicateVertices(currentMesh);
    displayMesh(currentMesh);
    console.log('Removed', removed, 'duplicate vertices');
  });

  document.getElementById('btn-repair-mesh').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    geomlib.repairMesh(currentMesh, false); // false = don't fill holes
    displayMesh(currentMesh);
    console.log('Mesh repaired');

    // Re-analyze to show improvements
    const issues = geomlib.analyzeMeshTopology(currentMesh);
    displayHealingStats(issues);
  });

  // Remeshing
  document.getElementById('btn-remesh-uniform').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    const targetLength = geomlib.computeMeanEdgeLength(currentMesh);
    geomlib.remeshUniform(currentMesh, targetLength, 10);
    displayMesh(currentMesh);

    const stats = geomlib.getRemeshingStats(currentMesh);
    displayRemeshingStats(stats);
    console.log('Uniform remeshing complete:', stats);
  });

  document.getElementById('btn-remesh-adaptive').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    const meanLength = geomlib.computeMeanEdgeLength(currentMesh);
    geomlib.remeshAdaptive(currentMesh, meanLength * 0.5, meanLength * 2.0, 10);
    displayMesh(currentMesh);

    const stats = geomlib.getRemeshingStats(currentMesh);
    displayRemeshingStats(stats);
    console.log('Adaptive remeshing complete:', stats);
  });

  document.getElementById('btn-remesh-custom').addEventListener('click', () => {
    if (!currentMesh) {
      alert('Create a mesh first!');
      return;
    }

    const targetLength = parseFloat(document.getElementById('remesh-length-slider').value);
    geomlib.remeshUniform(currentMesh, targetLength, 10);
    displayMesh(currentMesh);

    const stats = geomlib.getRemeshingStats(currentMesh);
    displayRemeshingStats(stats);
    console.log('Custom remeshing complete:', stats);
  });

  document.getElementById('remesh-length-slider').addEventListener('input', (e) => {
    document.getElementById('remesh-length-value').textContent = e.target.value;
  });

  // Display options
  document.getElementById('chk-wireframe').addEventListener('change', (e) => {
    state.showWireframe = e.target.checked;
    if (currentMesh) displayMesh(currentMesh);
    if (currentProxy) displayMesh(currentProxy, true);
  });

  document.getElementById('chk-proxy').addEventListener('change', (e) => {
    state.showProxy = e.target.checked;
    updateProxyVisibility();
  });

  document.getElementById('chk-bvh').addEventListener('change', (e) => {
    state.showBVH = e.target.checked;

    // Show/hide depth control
    const depthControl = document.getElementById('bvh-depth-control');
    depthControl.style.display = state.showBVH ? 'block' : 'none';

    visualizeBVH();
  });

  document.getElementById('bvh-depth-slider').addEventListener('input', (e) => {
    state.bvhMaxDepth = parseInt(e.target.value);
    document.getElementById('bvh-depth-value').textContent = state.bvhMaxDepth;
    visualizeBVH();
  });

  document.getElementById('chk-points').addEventListener('change', (e) => {
    state.showPoints = e.target.checked;
    updatePointCloudVisibility();
  });

  // File upload handler
  document.getElementById('file-input').addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;

    try {
      console.log(`Loading file: ${file.name}`);
      const result = await loadPointCloudFile(file);

      if (result.type === 'ply') {
        // Create point cloud from PLY data
        const pc = new geomlib.PointCloud();
        const data = result.data;

        for (let i = 0; i < data.vertexCount; i++) {
          const x = data.positions[i * 3 + 0];
          const y = data.positions[i * 3 + 1];
          const z = data.positions[i * 3 + 2];

          if (data.normals) {
            const nx = data.normals[i * 3 + 0];
            const ny = data.normals[i * 3 + 1];
            const nz = data.normals[i * 3 + 2];
            pc.addPointWithNormal(x, y, z, nx, ny, nz);
          } else {
            pc.addPoint(x, y, z);
          }
        }

        currentPointCloud = pc;
        visualizePointCloud(currentPointCloud);
        state.showPoints = true;
        document.getElementById('chk-points').checked = true;
        updatePointCloudVisibility();
        clearMesh();

        console.log(`Loaded ${data.vertexCount} points from PLY file`);
      } else if (result.type === 'splat') {
        // Create point cloud from Gaussian splat data
        const data = result.data;
        const positions = Array.from(data.positions);

        currentPointCloud = geomlib.gaussiansSamplePoints(positions, 3);
        visualizePointCloud(currentPointCloud);
        state.showPoints = true;
        document.getElementById('chk-points').checked = true;
        updatePointCloudVisibility();
        clearMesh();

        console.log(`Loaded ${data.splatCount} Gaussian splats, sampled to ${currentPointCloud.getPointCount()} points`);
      }
    } catch (error) {
      console.error('Error loading file:', error);
      alert(`Error loading file: ${error.message}`);
    }

    // Clear file input
    e.target.value = '';
  });
}

function clearMesh() {
  if (mainMeshGroup) {
    scene.remove(mainMeshGroup);
    mainMeshGroup = null;
  }
  currentMesh = null;
  document.getElementById('mesh-stats').innerHTML = '';
}

function clearProxy() {
  if (proxyMeshGroup) {
    scene.remove(proxyMeshGroup);
    proxyMeshGroup = null;
  }
  currentProxy = null;
  state.showProxy = false;
  document.getElementById('chk-proxy').checked = false;
}

function clearBVH() {
  if (currentBVH) {
    currentBVH = null;
    document.getElementById('bvh-stats').innerHTML = '';
  }
  clearRayLines();
  clearBVHBoxes();
  state.showBVH = false;
  document.getElementById('chk-bvh').checked = false;
  document.getElementById('chk-bvh').disabled = true;
  document.getElementById('bvh-depth-control').style.display = 'none';
}

function updateProxyVisibility() {
  if (proxyMeshGroup) {
    proxyMeshGroup.visible = state.showProxy;
  }
  document.getElementById('chk-proxy').checked = state.showProxy;
}

function clearRayLines() {
  if (rayLines) {
    scene.remove(rayLines);
    rayLines = null;
  }
}

// ============================================================================
// BVH Visualization
// ============================================================================

function visualizeBVH() {
  clearBVHBoxes();

  if (!currentBVH || !state.showBVH) {
    return;
  }

  const boxes = currentBVH.getBoundingBoxes(state.bvhMaxDepth);
  console.log(`Visualizing ${boxes.length} BVH nodes at depth <= ${state.bvhMaxDepth}`);

  bvhBoxesGroup = new THREE.Group();

  // Color map for different depths
  const depthColors = [
    0xff0000, // Red - root
    0xff8800, // Orange
    0xffff00, // Yellow
    0x00ff00, // Green
    0x00ffff, // Cyan
    0x0088ff, // Light Blue
    0x0000ff, // Blue
    0x8800ff, // Purple
    0xff00ff, // Magenta
    0xff0088, // Pink
  ];

  for (let i = 0; i < boxes.length; i++) {
    const box = boxes[i];

    const width = box.maxX - box.minX;
    const height = box.maxY - box.minY;
    const depth = box.maxZ - box.minZ;

    const centerX = (box.minX + box.maxX) / 2;
    const centerY = (box.minY + box.maxY) / 2;
    const centerZ = (box.minZ + box.maxZ) / 2;

    // Create box geometry
    const geometry = new THREE.BoxGeometry(width, height, depth);

    // Color based on depth
    const color = depthColors[box.depth % depthColors.length];
    const material = new THREE.MeshBasicMaterial({
      color: color,
      wireframe: true,
      transparent: true,
      opacity: box.isLeaf ? 0.6 : 0.3
    });

    const boxMesh = new THREE.Mesh(geometry, material);
    boxMesh.position.set(centerX, centerY, centerZ);

    bvhBoxesGroup.add(boxMesh);
  }

  scene.add(bvhBoxesGroup);
}

function clearBVHBoxes() {
  if (bvhBoxesGroup) {
    scene.remove(bvhBoxesGroup);
    bvhBoxesGroup = null;
  }
}

// ============================================================================
// Point Cloud Visualization
// ============================================================================

function visualizePointCloud(pointCloud) {
  clearPointCloud();

  const points = pointCloud.getPoints();
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(points, 3));

  const material = new THREE.PointsMaterial({
    color: 0x00ffff,
    size: 0.05,
    sizeAttenuation: true
  });

  pointCloudPoints = new THREE.Points(geometry, material);
  pointCloudPoints.visible = state.showPoints;
  scene.add(pointCloudPoints);

  // Update stats
  const stats = pointCloud.getStats();
  const statsDiv = document.getElementById('pointcloud-stats');
  statsDiv.innerHTML = `
    <strong>Point Cloud:</strong><br>
    Points: ${stats.pointCount}<br>
    Avg Spacing: ${stats.avgSpacing.toFixed(4)}<br>
    Density: ${stats.density.toFixed(2)}
  `;
}

function clearPointCloud() {
  if (pointCloudPoints) {
    scene.remove(pointCloudPoints);
    pointCloudPoints = null;
  }
  document.getElementById('pointcloud-stats').innerHTML = '';
}

function updatePointCloudVisibility() {
  if (pointCloudPoints) {
    pointCloudPoints.visible = state.showPoints;
  }
}

// ============================================================================
// Ray Tracing Visualization
// ============================================================================

function castRandomRays(count) {
  clearRayLines();

  const positions = [];
  const colors = [];

  for (let i = 0; i < count; i++) {
    // Random ray origin on sphere around mesh
    const theta = Math.random() * Math.PI * 2;
    const phi = Math.random() * Math.PI;
    const radius = 3;

    const ox = radius * Math.sin(phi) * Math.cos(theta);
    const oy = radius * Math.sin(phi) * Math.sin(theta);
    const oz = radius * Math.cos(phi);

    // Ray direction towards origin
    const dx = -ox;
    const dy = -oy;
    const dz = -oz;
    const len = Math.sqrt(dx * dx + dy * dy + dz * dz);

    // Trace ray
    const result = currentBVH.traceRay(ox, oy, oz, dx / len, dy / len, dz / len);

    // Add line
    positions.push(ox, oy, oz);

    if (result.hit) {
      positions.push(result.x, result.y, result.z);
      colors.push(0, 1, 0, 0, 1, 0); // Green for hit
    } else {
      positions.push(ox + dx, oy + dy, oz + dz);
      colors.push(1, 0, 0, 1, 0, 0); // Red for miss
    }
  }

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

  const material = new THREE.LineBasicMaterial({
    vertexColors: true,
    transparent: true,
    opacity: 0.6
  });

  rayLines = new THREE.LineSegments(geometry, material);
  scene.add(rayLines);
}

// ============================================================================
// Initialization
// ============================================================================

async function init() {
  console.log('Initializing Geomlib WebAssembly demo...');

  // Initialize Three.js
  initThreeJS();

  // Load WebAssembly module
  try {
    geomlib = await createGeomlibModule();
    console.log('Geomlib module loaded successfully');

    // Setup UI event handlers
    setupEventHandlers();

    // Hide loading screen
    const loading = document.getElementById('loading');
    loading.classList.add('hidden');

    // Create initial sphere
    currentMesh = geomlib.createSphere(1.0, 32, 16);
    displayMesh(currentMesh);

    console.log('Demo ready!');
  } catch (error) {
    console.error('Failed to load Geomlib module:', error);
    alert('Failed to load WebAssembly module. Check console for details.');
  }
}

// Start the application
init();
