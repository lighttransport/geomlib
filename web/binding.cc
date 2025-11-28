/*
 * geomlib - WebAssembly/JavaScript Bindings
 * Using Emscripten embind
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cstdlib>

extern "C" {
    #include "geomlib/geomlib.h"
}

using namespace emscripten;

// ============================================================================
// Helper Classes for JavaScript API
// ============================================================================

class Vec3 {
public:
    float x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    geomlib_vec3 to_c() const {
        return geomlib_vec3_create(x, y, z);
    }

    static Vec3 from_c(geomlib_vec3 v) {
        return Vec3(v.x, v.y, v.z);
    }
};

class Mesh {
public:
    geomlib_trimesh mesh;

    Mesh() {
        geomlib_trimesh_init(&mesh, 64, 64);
    }

    ~Mesh() {
        geomlib_trimesh_free(&mesh);
    }

    // Prevent copying
    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    void addVertex(float x, float y, float z) {
        geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(x, y, z), nullptr, nullptr);
    }

    void addTriangle(uint32_t i0, uint32_t i1, uint32_t i2) {
        geomlib_trimesh_add_triangle(&mesh, i0, i1, i2);
    }

    int getVertexCount() const {
        return mesh.vertex_count;
    }

    int getTriangleCount() const {
        return mesh.triangle_count;
    }

    // Get vertices as flat array for Three.js
    val getVertices() const {
        val vertices = val::array();
        for (size_t i = 0; i < mesh.vertex_count; i++) {
            vertices.call<void>("push", mesh.vertices[i].x);
            vertices.call<void>("push", mesh.vertices[i].y);
            vertices.call<void>("push", mesh.vertices[i].z);
        }
        return vertices;
    }

    // Get indices as flat array for Three.js
    val getIndices() const {
        val indices = val::array();
        for (size_t i = 0; i < mesh.triangle_count * 3; i++) {
            indices.call<void>("push", mesh.indices[i]);
        }
        return indices;
    }

    // Get normals (compute if needed)
    val getNormals() const {
        val normals = val::array();

        // Compute per-vertex normals if not present
        if (!mesh.normals) {
            geomlib_trimesh_compute_normals(const_cast<geomlib_trimesh*>(&mesh));
        }

        if (mesh.normals) {
            for (size_t i = 0; i < mesh.vertex_count; i++) {
                normals.call<void>("push", mesh.normals[i].x);
                normals.call<void>("push", mesh.normals[i].y);
                normals.call<void>("push", mesh.normals[i].z);
            }
        }

        return normals;
    }

    float getSurfaceArea() const {
        return geomlib_trimesh_surface_area(&mesh);
    }

    float getVolume() const {
        return geomlib_trimesh_volume(&mesh);
    }

    Vec3 getCentroid() const {
        return Vec3::from_c(geomlib_trimesh_centroid(&mesh));
    }
};

// ============================================================================
// Mesh Creation Functions
// ============================================================================

Mesh* createBox(float width, float height, float depth) {
    Mesh* m = new Mesh();
    geomlib_create_box(&m->mesh, width, height, depth);
    return m;
}

Mesh* createSphere(float radius, int segments, int rings) {
    Mesh* m = new Mesh();
    geomlib_create_sphere(&m->mesh, radius, segments, rings);
    return m;
}

Mesh* createCylinder(float radius, float height, int segments) {
    Mesh* m = new Mesh();
    geomlib_create_cylinder(&m->mesh, radius, height, segments);
    return m;
}

Mesh* createPlane(float width, float height, int subdivisions) {
    Mesh* m = new Mesh();
    geomlib_create_plane(&m->mesh, width, height, subdivisions);
    return m;
}

// ============================================================================
// Collision Proxy Generation
// ============================================================================

Mesh* generateCollisionProxy(const Mesh& mesh, int strategy, int targetTriangles) {
    Mesh* proxy = new Mesh();

    geomlib_proxy_strategy strat;
    switch (strategy) {
        case 0: strat = GEOMLIB_PROXY_CONVEX_HULL; break;
        case 1: strat = GEOMLIB_PROXY_SIMPLIFIED; break;
        case 2: strat = GEOMLIB_PROXY_OBB; break;
        case 3: strat = GEOMLIB_PROXY_AABB; break;
        default: strat = GEOMLIB_PROXY_SIMPLIFIED; break;
    }

    geomlib_generate_collision_proxy(&mesh.mesh, strat, targetTriangles, &proxy->mesh);
    return proxy;
}

Mesh* simplifyMesh(const Mesh& mesh, float reductionRatio) {
    Mesh* simplified = new Mesh();

    // Copy mesh
    simplified->mesh.vertex_count = mesh.mesh.vertex_count;
    simplified->mesh.triangle_count = mesh.mesh.triangle_count;
    simplified->mesh.vertex_capacity = mesh.mesh.vertex_count;
    simplified->mesh.triangle_capacity = mesh.mesh.triangle_count;

    simplified->mesh.vertices = (geomlib_vec3*)malloc(mesh.mesh.vertex_count * sizeof(geomlib_vec3));
    simplified->mesh.indices = (uint32_t*)malloc(mesh.mesh.triangle_count * 3 * sizeof(uint32_t));

    memcpy(simplified->mesh.vertices, mesh.mesh.vertices, mesh.mesh.vertex_count * sizeof(geomlib_vec3));
    memcpy(simplified->mesh.indices, mesh.mesh.indices, mesh.mesh.triangle_count * 3 * sizeof(uint32_t));

    simplified->mesh.normals = nullptr;
    simplified->mesh.texcoords = nullptr;

    geomlib_simplify_mesh_ratio(&simplified->mesh, reductionRatio, 0);
    return simplified;
}

Mesh* computeConvexHull(const Mesh& mesh) {
    Mesh* hull = new Mesh();
    geomlib_compute_convex_hull(&mesh.mesh, &hull->mesh);
    return hull;
}

// ============================================================================
// BVH and Ray Tracing
// ============================================================================

class BVH {
public:
    geomlib_bvh bvh;
    const Mesh* mesh;

    BVH(const Mesh* m) : mesh(m) {
        geomlib_bvh_build(&bvh, &m->mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    }

    ~BVH() {
        geomlib_bvh_free(&bvh);
    }

    // Prevent copying
    BVH(const BVH&) = delete;
    BVH& operator=(const BVH&) = delete;

    val traceRay(float ox, float oy, float oz, float dx, float dy, float dz) {
        geomlib_ray ray;
        ray.origin = geomlib_vec3_create(ox, oy, oz);
        ray.direction = geomlib_vec3_create(dx, dy, dz);
        ray.t_min = 0.0f;
        ray.t_max = 1000.0f;

        geomlib_ray_hit hit;
        bool intersects = geomlib_bvh_trace_ray(&bvh, &ray, &hit);

        val result = val::object();
        result.set("hit", intersects && hit.hit);

        if (intersects && hit.hit) {
            result.set("t", hit.t);
            result.set("x", hit.position.x);
            result.set("y", hit.position.y);
            result.set("z", hit.position.z);
            result.set("nx", hit.normal.x);
            result.set("ny", hit.normal.y);
            result.set("nz", hit.normal.z);
        }

        return result;
    }

    val getStats() {
        uint32_t max_depth, leaf_count;
        float avg_leaf_size;
        geomlib_bvh_get_stats(&bvh, &max_depth, &leaf_count, &avg_leaf_size);

        val stats = val::object();
        stats.set("maxDepth", max_depth);
        stats.set("leafCount", leaf_count);
        stats.set("avgLeafSize", avg_leaf_size);
        return stats;
    }

    val getBoundingBoxes(int maxDepth = -1) {
        val boxes = val::array();

        if (bvh.node_count == 0) return boxes;

        // Traverse BVH and collect bounding boxes
        collectBoxes(boxes, 0, 0, maxDepth);

        return boxes;
    }

private:
    void collectBoxes(val& boxes, size_t nodeIndex, int depth, int maxDepth) {
        if (nodeIndex >= bvh.node_count) return;
        if (maxDepth >= 0 && depth > maxDepth) return;

        const geomlib_bvh_node* node = &bvh.nodes[nodeIndex];

        // Add this node's bounding box
        val box = val::object();
        box.set("minX", node->bounds.min.x);
        box.set("minY", node->bounds.min.y);
        box.set("minZ", node->bounds.min.z);
        box.set("maxX", node->bounds.max.x);
        box.set("maxY", node->bounds.max.y);
        box.set("maxZ", node->bounds.max.z);
        box.set("isLeaf", node->is_leaf);
        box.set("depth", depth);

        boxes.call<void>("push", box);

        // Recurse to children if interior node
        if (!node->is_leaf) {
            collectBoxes(boxes, node->data.interior.left_child, depth + 1, maxDepth);
            collectBoxes(boxes, node->data.interior.right_child, depth + 1, maxDepth);
        }
    }

public:
};

BVH* buildBVH(const Mesh& mesh) {
    return new BVH(&mesh);
}

// ============================================================================
// Mesh Analysis
// ============================================================================

val analyzeMeshQuality(const Mesh& mesh) {
    geomlib_mesh_quality quality;
    geomlib_analyze_mesh_quality(&mesh.mesh, &quality);

    val result = val::object();
    result.set("minEdgeLength", quality.min_edge_length);
    result.set("maxEdgeLength", quality.max_edge_length);
    result.set("avgEdgeLength", quality.avg_edge_length);
    result.set("minTriangleArea", quality.min_triangle_area);
    result.set("maxTriangleArea", quality.max_triangle_area);
    result.set("degenerateCount", (int)quality.degenerate_count);

    return result;
}

float estimateCollisionCost(const Mesh& mesh) {
    return geomlib_estimate_collision_cost(&mesh.mesh);
}

// ============================================================================
// Point Cloud and Surface Reconstruction
// ============================================================================

class PointCloud {
public:
    geomlib_pointcloud pc;

    PointCloud() {
        geomlib_pointcloud_init(&pc, 1024);
    }

    ~PointCloud() {
        geomlib_pointcloud_free(&pc);
    }

    // Prevent copying
    PointCloud(const PointCloud&) = delete;
    PointCloud& operator=(const PointCloud&) = delete;

    void addPoint(float x, float y, float z) {
        geomlib_vec3 pos = geomlib_vec3_create(x, y, z);
        geomlib_pointcloud_add_point(&pc, pos, nullptr, 1.0f);
    }

    void addPointWithNormal(float x, float y, float z, float nx, float ny, float nz) {
        geomlib_vec3 pos = geomlib_vec3_create(x, y, z);
        geomlib_vec3 normal = geomlib_vec3_create(nx, ny, nz);
        geomlib_pointcloud_add_point(&pc, pos, &normal, 1.0f);
    }

    int getPointCount() const {
        return (int)pc.point_count;
    }

    val getPoints() const {
        val points = val::array();
        for (size_t i = 0; i < pc.point_count; i++) {
            points.call<void>("push", pc.points[i].x);
            points.call<void>("push", pc.points[i].y);
            points.call<void>("push", pc.points[i].z);
        }
        return points;
    }

    val getStats() const {
        float avg_spacing, density;
        geomlib_pointcloud_get_stats(&pc, &avg_spacing, &density);

        val stats = val::object();
        stats.set("pointCount", (int)pc.point_count);
        stats.set("avgSpacing", avg_spacing);
        stats.set("density", density);
        return stats;
    }
};

// Create point cloud from Gaussian splat positions
PointCloud* gaussiansToPointCloud(val positions_js) {
    std::vector<float> positions = convertJSArrayToNumberVector<float>(positions_js);

    if (positions.size() % 3 != 0) {
        return nullptr;
    }

    PointCloud* pc = new PointCloud();
    geomlib_gaussians_to_pointcloud(positions.data(), positions.size() / 3, &pc->pc);

    return pc;
}

// Sample points from Gaussians
PointCloud* gaussiansSamplePoints(val positions_js, int samplesPerGaussian) {
    std::vector<float> positions = convertJSArrayToNumberVector<float>(positions_js);

    if (positions.size() % 3 != 0) {
        return nullptr;
    }

    PointCloud* pc = new PointCloud();
    geomlib_gaussians_sample_points(
        positions.data(), nullptr, nullptr,
        positions.size() / 3, samplesPerGaussian, &pc->pc
    );

    return pc;
}

// Reconstruct surface from point cloud
Mesh* reconstructSurface(const PointCloud& pc, int method, float voxelSize) {
    Mesh* mesh = new Mesh();

    geomlib_surface_params params;
    params.method = (geomlib_surface_method)method;
    params.voxel_size = voxelSize;
    params.max_edge_length = 0.0f;  // Auto
    params.max_neighbors = 20;
    params.search_radius = 0.0f;  // Auto

    geomlib_reconstruct_surface(&pc.pc, &params, &mesh->mesh);

    return mesh;
}

// Quick auto surface reconstruction
Mesh* reconstructSurfaceAuto(const PointCloud& pc, int targetFaces) {
    Mesh* mesh = new Mesh();
    geomlib_reconstruct_surface_auto(&pc.pc, targetFaces, &mesh->mesh);
    return mesh;
}

// ============================================================================
// Mesh Healing
// ============================================================================

void fixFaceOrientation(Mesh& mesh, int strategy) {
    geomlib_orientation_params params = geomlib_orientation_params_default(
        (geomlib_orientation_strategy)strategy
    );
    geomlib_fix_face_orientation(&mesh.mesh, &params, nullptr);
}

int removeDegenerateTriangles(Mesh& mesh) {
    return geomlib_remove_degenerate_triangles(&mesh.mesh, 1e-10f);
}

int removeDuplicateVertices(Mesh& mesh) {
    return geomlib_remove_duplicate_vertices(&mesh.mesh, 1e-6f);
}

void repairMesh(Mesh& mesh, bool fillHoles) {
    geomlib_repair_mesh(&mesh.mesh, fillHoles);
}

val analyzeMeshTopology(const Mesh& mesh) {
    geomlib_mesh_issues issues;
    geomlib_analyze_mesh_topology(&mesh.mesh, &issues);

    val result = val::object();
    result.set("degenerateTriangles", (int)issues.degenerate_triangles);
    result.set("duplicateVertices", (int)issues.duplicate_vertices);
    result.set("boundaryEdges", (int)issues.boundary_edges);
    result.set("nonManifoldEdges", (int)issues.non_manifold_edges);
    result.set("nonManifoldVertices", (int)issues.non_manifold_vertices);
    result.set("inconsistentFaces", (int)issues.inconsistent_faces);
    result.set("isManifold", issues.is_manifold);
    result.set("isWatertight", issues.is_watertight);
    result.set("isOriented", issues.is_oriented);

    return result;
}

// Remeshing functions
void remeshUniform(Mesh& mesh, float targetEdgeLength, int iterations) {
    geomlib_remesh_uniform(&mesh.mesh, targetEdgeLength, iterations);
}

void remeshAdaptive(Mesh& mesh, float minEdgeLength, float maxEdgeLength, int iterations) {
    geomlib_remesh_adaptive(&mesh.mesh, minEdgeLength, maxEdgeLength, iterations);
}

float computeTargetEdgeLength(const Mesh& mesh, int targetVertexCount) {
    return geomlib_compute_target_edge_length(&mesh.mesh, targetVertexCount);
}

float computeMeanEdgeLength(const Mesh& mesh) {
    return geomlib_compute_mean_edge_length(&mesh.mesh);
}

val getRemeshingStats(const Mesh& mesh) {
    geomlib_remeshing_stats stats;
    geomlib_compute_remeshing_stats(&mesh.mesh, &stats);

    val result = val::object();
    result.set("edgesSplit", stats.edges_split);
    result.set("edgesCollapsed", stats.edges_collapsed);
    result.set("edgesFlipped", stats.edges_flipped);
    result.set("minEdgeLength", stats.min_edge_length);
    result.set("maxEdgeLength", stats.max_edge_length);
    result.set("meanEdgeLength", stats.mean_edge_length);
    result.set("minAngle", stats.min_angle);
    result.set("maxAngle", stats.max_angle);
    result.set("irregularVertices", stats.irregular_vertices);

    return result;
}

// Parameterization functions
val parameterizeHarmonic(const Mesh& mesh, int boundaryShape) {
    geomlib_uv_coords uvs = geomlib_uv_coords_create(mesh.mesh.vertex_count, mesh.mesh.triangle_count);

    geomlib_boundary_shape shape = (geomlib_boundary_shape)boundaryShape;
    int result = geomlib_parameterize_harmonic(&mesh.mesh, shape, &uvs);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val uvArray = val::array();
        for (size_t i = 0; i < uvs.vertex_count; i++) {
            uvArray.call<void>("push", uvs.uvs[i].x);
            uvArray.call<void>("push", uvs.uvs[i].y);
        }
        output.set("uvs", uvArray);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    geomlib_uv_coords_free(&uvs);
    return output;
}

val parameterizeLSCM(const Mesh& mesh, int boundaryShape) {
    geomlib_uv_coords uvs = geomlib_uv_coords_create(mesh.mesh.vertex_count, mesh.mesh.triangle_count);

    geomlib_boundary_shape shape = (geomlib_boundary_shape)boundaryShape;
    int result = geomlib_parameterize_lscm(&mesh.mesh, shape, &uvs);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val uvArray = val::array();
        for (size_t i = 0; i < uvs.vertex_count; i++) {
            uvArray.call<void>("push", uvs.uvs[i].x);
            uvArray.call<void>("push", uvs.uvs[i].y);
        }
        output.set("uvs", uvArray);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    geomlib_uv_coords_free(&uvs);
    return output;
}

// ============================================================================
// Curvature Computation
// ============================================================================

val computeCurvature(const Mesh& mesh) {
    geomlib_mesh_curvature curvature = geomlib_curvature_create(mesh.mesh.vertex_count);

    int result = geomlib_compute_curvature(&mesh.mesh, &curvature);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val meanCurvatures = val::array();
        val gaussianCurvatures = val::array();
        val principalMax = val::array();
        val principalMin = val::array();

        for (size_t i = 0; i < curvature.vertex_count; i++) {
            meanCurvatures.call<void>("push", curvature.curvatures[i].mean);
            gaussianCurvatures.call<void>("push", curvature.curvatures[i].gaussian);
            principalMax.call<void>("push", curvature.curvatures[i].principal_max);
            principalMin.call<void>("push", curvature.curvatures[i].principal_min);
        }

        output.set("mean", meanCurvatures);
        output.set("gaussian", gaussianCurvatures);
        output.set("principalMax", principalMax);
        output.set("principalMin", principalMin);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    geomlib_curvature_free(&curvature);
    return output;
}

val computeMeanCurvature(const Mesh& mesh) {
    float *mean_curvatures = (float*)malloc(mesh.mesh.vertex_count * sizeof(float));

    int result = geomlib_compute_mean_curvature(&mesh.mesh, mean_curvatures);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val curvatures = val::array();
        for (size_t i = 0; i < mesh.mesh.vertex_count; i++) {
            curvatures.call<void>("push", mean_curvatures[i]);
        }
        output.set("curvatures", curvatures);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    free(mean_curvatures);
    return output;
}

val computeGaussianCurvature(const Mesh& mesh) {
    float *gaussian_curvatures = (float*)malloc(mesh.mesh.vertex_count * sizeof(float));

    int result = geomlib_compute_gaussian_curvature(&mesh.mesh, gaussian_curvatures);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val curvatures = val::array();
        for (size_t i = 0; i < mesh.mesh.vertex_count; i++) {
            curvatures.call<void>("push", gaussian_curvatures[i]);
        }
        output.set("curvatures", curvatures);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    free(gaussian_curvatures);
    return output;
}

// ============================================================================
// Geodesic Distances
// ============================================================================

val computeGeodesicDistances(const Mesh& mesh, uint32_t source_vertex) {
    geomlib_geodesic_distances geo = geomlib_geodesic_create(mesh.mesh.vertex_count);

    int result = geomlib_compute_geodesic_from_vertex(&mesh.mesh, source_vertex, &geo);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val distances = val::array();
        for (size_t i = 0; i < geo.vertex_count; i++) {
            distances.call<void>("push", geo.distances[i]);
        }
        output.set("distances", distances);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    geomlib_geodesic_free(&geo);
    return output;
}

val computeGeodesicDistancesMultiple(const Mesh& mesh, val source_vertices_js) {
    // Convert JavaScript array to C array
    unsigned int length = source_vertices_js["length"].as<unsigned int>();
    uint32_t *source_vertices = (uint32_t*)malloc(length * sizeof(uint32_t));

    for (unsigned int i = 0; i < length; i++) {
        source_vertices[i] = source_vertices_js[i].as<uint32_t>();
    }

    geomlib_geodesic_distances geo = geomlib_geodesic_create(mesh.mesh.vertex_count);

    int result = geomlib_compute_geodesic_heat(&mesh.mesh, source_vertices, length, 0.0f, &geo);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        val distances = val::array();
        for (size_t i = 0; i < geo.vertex_count; i++) {
            distances.call<void>("push", geo.distances[i]);
        }
        output.set("distances", distances);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    free(source_vertices);
    geomlib_geodesic_free(&geo);
    return output;
}

// ============================================================================
// Mesh Decimation
// ============================================================================

val decimateQEM(Mesh& mesh, float target_ratio) {
    geomlib_decimation_stats stats;
    geomlib_decimation_params params = geomlib_decimation_params_default(1.0f - target_ratio);

    int result = geomlib_decimate_qem(&mesh.mesh, &params, &stats);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        output.set("originalVertices", (unsigned int)stats.original_vertices);
        output.set("originalTriangles", (unsigned int)stats.original_triangles);
        output.set("finalVertices", (unsigned int)stats.final_vertices);
        output.set("finalTriangles", (unsigned int)stats.final_triangles);
        output.set("edgesCollapsed", (unsigned int)stats.edges_collapsed);
        output.set("reductionRatio", stats.reduction_ratio);
        output.set("maxError", stats.max_error);
        output.set("meanError", stats.mean_error);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    return output;
}

val decimateToVertexCount(Mesh& mesh, uint32_t target_vertices) {
    geomlib_decimation_stats stats;
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_vertex_count = target_vertices;

    int result = geomlib_decimate_qem(&mesh.mesh, &params, &stats);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        output.set("originalVertices", (unsigned int)stats.original_vertices);
        output.set("finalVertices", (unsigned int)stats.final_vertices);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    return output;
}

// ============================================================================
// Mesh Smoothing
// ============================================================================

val smoothLaplacian(Mesh& mesh, float lambda, int iterations) {
    int result = geomlib_smooth_laplacian(&mesh.mesh, lambda, iterations);

    val output = val::object();
    output.set("success", result == GEOMLIB_MESH_SUCCESS);
    return output;
}

val smoothTaubin(Mesh& mesh, float lambda, float mu, int iterations) {
    int result = geomlib_smooth_taubin(&mesh.mesh, lambda, mu, iterations);

    val output = val::object();
    output.set("success", result == GEOMLIB_MESH_SUCCESS);
    return output;
}

val smoothNormals(Mesh& mesh, int iterations) {
    int result = geomlib_smooth_normals(&mesh.mesh, iterations);

    val output = val::object();
    output.set("success", result == GEOMLIB_MESH_SUCCESS);
    return output;
}

// ============================================================================
// Mesh Quality Analysis
// ============================================================================

val computeMeshQuality(const Mesh& mesh) {
    geomlib_mesh_quality_stats stats;
    geomlib_quality_thresholds thresholds = geomlib_quality_thresholds_default();

    int result = geomlib_compute_mesh_quality(&mesh.mesh, &thresholds, &stats);

    val output = val::object();
    if (result == GEOMLIB_MESH_SUCCESS) {
        output.set("minTriangleArea", stats.min_triangle_area);
        output.set("maxTriangleArea", stats.max_triangle_area);
        output.set("meanTriangleArea", stats.mean_triangle_area);

        output.set("minAngle", stats.min_angle_global);
        output.set("maxAngle", stats.max_angle_global);
        output.set("meanMinAngle", stats.mean_min_angle);

        output.set("minAspectRatio", stats.min_aspect_ratio);
        output.set("maxAspectRatio", stats.max_aspect_ratio);
        output.set("meanAspectRatio", stats.mean_aspect_ratio);

        output.set("minShapeMeasure", stats.min_shape_measure);
        output.set("meanShapeMeasure", stats.mean_shape_measure);

        output.set("minEdgeLength", stats.min_edge_length);
        output.set("maxEdgeLength", stats.max_edge_length);
        output.set("meanEdgeLength", stats.mean_edge_length);

        output.set("degenerateCount", (unsigned int)stats.degenerate_count);
        output.set("sliverCount", (unsigned int)stats.sliver_count);
        output.set("boundaryEdges", (unsigned int)stats.boundary_edges);

        output.set("minValence", stats.min_valence);
        output.set("maxValence", stats.max_valence);
        output.set("meanValence", stats.mean_valence);
        output.set("irregularVertices", (unsigned int)stats.irregular_vertices);

        output.set("qualityScore", stats.quality_score);
        output.set("success", true);
    } else {
        output.set("success", false);
    }

    return output;
}

val checkMeshQuality(const Mesh& mesh) {
    geomlib_quality_thresholds thresholds = geomlib_quality_thresholds_default();
    bool passed = geomlib_check_mesh_quality(&mesh.mesh, &thresholds);

    val output = val::object();
    output.set("passed", passed);
    return output;
}

// ============================================================================
// Texture Atlas Generation
// ============================================================================

val generateTextureAtlas(Mesh& mesh) {
    geomlib_atlas_params params = geomlib_atlas_params_default();
    geomlib_atlas atlas = geomlib_atlas_create(256);

    int result = geomlib_generate_atlas(&mesh.mesh, &params, &atlas);

    val output = val::object();

    if (result == GEOMLIB_MESH_SUCCESS) {
        output.set("success", true);
        output.set("chartCount", (unsigned int)atlas.chart_count);
        output.set("coverage", atlas.coverage);
        output.set("avgStretch", atlas.avg_stretch);
        output.set("maxStretch", atlas.max_stretch);

        // Return chart info
        val charts = val::array();
        for (size_t c = 0; c < atlas.chart_count; c++) {
            val chart = val::object();
            chart.set("triangleCount", (unsigned int)atlas.charts[c].triangle_count);
            chart.set("vertexCount", (unsigned int)atlas.charts[c].vertex_count);

            val offset = val::object();
            offset.set("x", atlas.charts[c].atlas_offset.x);
            offset.set("y", atlas.charts[c].atlas_offset.y);
            chart.set("offset", offset);

            val bounds = val::object();
            bounds.set("minX", atlas.charts[c].uv_min.x);
            bounds.set("minY", atlas.charts[c].uv_min.y);
            bounds.set("maxX", atlas.charts[c].uv_max.x);
            bounds.set("maxY", atlas.charts[c].uv_max.y);
            chart.set("bounds", bounds);

            charts.call<void>("push", chart);
        }
        output.set("charts", charts);

        // Return UVs array
        if (mesh.mesh.texcoords) {
            val uvs = val::array();
            for (size_t i = 0; i < mesh.mesh.vertex_count; i++) {
                uvs.call<void>("push", mesh.mesh.texcoords[i].x);
                uvs.call<void>("push", mesh.mesh.texcoords[i].y);
            }
            output.set("uvs", uvs);
        }
    } else {
        output.set("success", false);
    }

    geomlib_atlas_free(&atlas);
    return output;
}

val generateTextureAtlasCustom(
    Mesh& mesh,
    int segmentation,
    int packing,
    int parameterization,
    float gutterWidth
) {
    geomlib_atlas_params params = geomlib_atlas_params_default();
    params.segmentation = (geomlib_atlas_segmentation)segmentation;
    params.packing = (geomlib_atlas_packing)packing;
    params.parameterization = (geomlib_param_method)parameterization;
    params.gutter_width = gutterWidth;

    geomlib_atlas atlas = geomlib_atlas_create(256);

    int result = geomlib_generate_atlas(&mesh.mesh, &params, &atlas);

    val output = val::object();

    if (result == GEOMLIB_MESH_SUCCESS) {
        output.set("success", true);
        output.set("chartCount", (unsigned int)atlas.chart_count);
        output.set("coverage", atlas.coverage);

        // Return UVs array
        if (mesh.mesh.texcoords) {
            val uvs = val::array();
            for (size_t i = 0; i < mesh.mesh.vertex_count; i++) {
                uvs.call<void>("push", mesh.mesh.texcoords[i].x);
                uvs.call<void>("push", mesh.mesh.texcoords[i].y);
            }
            output.set("uvs", uvs);
        }
    } else {
        output.set("success", false);
    }

    geomlib_atlas_free(&atlas);
    return output;
}

// ============================================================================
// CFD Utilities
// ============================================================================

// Forward declaration
class TetraMesh;

class ScalarField {
public:
    geomlib_scalar_field field;

    ScalarField(int count) {
        field = geomlib_scalar_field_create(count);
    }

    ~ScalarField() {
        geomlib_scalar_field_free(&field);
    }

    ScalarField(const ScalarField&) = delete;
    ScalarField& operator=(const ScalarField&) = delete;

    void setValue(int index, float value) {
        if (index >= 0 && (size_t)index < field.count) {
            field.values[index] = value;
        }
    }

    float getValue(int index) const {
        if (index >= 0 && (size_t)index < field.count) {
            return field.values[index];
        }
        return 0.0f;
    }

    val getValues() const {
        val values = val::array();
        for (size_t i = 0; i < field.count; i++) {
            values.call<void>("push", field.values[i]);
        }
        return values;
    }
};

class VectorField {
public:
    geomlib_vector_field field;

    VectorField(int count) {
        field = geomlib_vector_field_create(count);
    }

    ~VectorField() {
        geomlib_vector_field_free(&field);
    }

    VectorField(const VectorField&) = delete;
    VectorField& operator=(const VectorField&) = delete;

    void setVector(int index, float x, float y, float z) {
        if (index >= 0 && (size_t)index < field.count) {
            field.values[index] = geomlib_vec3_create(x, y, z);
        }
    }

    val getVectors() const {
        val vectors = val::array();
        for (size_t i = 0; i < field.count; i++) {
            vectors.call<void>("push", field.values[i].x);
            vectors.call<void>("push", field.values[i].y);
            vectors.call<void>("push", field.values[i].z);
        }
        return vectors;
    }
};

VectorField* computeGradient(const Mesh& mesh, const ScalarField& field, int method) {
    VectorField* gradient = new VectorField(mesh.mesh.vertex_count);

    int result = geomlib_compute_gradient_scalar(
        &mesh.mesh,
        &field.field,
        (geomlib_gradient_method)method,
        &gradient->field
    );

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete gradient;
        return nullptr;
    }

    return gradient;
}

val traceStreamline(
    const Mesh& mesh,
    const VectorField& velocity,
    float seedX, float seedY, float seedZ,
    float stepSize, int maxSteps
) {
    geomlib_vec3 seed = geomlib_vec3_create(seedX, seedY, seedZ);

    geomlib_streamline_params params;
    params.method = GEOMLIB_INTEGRATE_RK4;
    params.step_size = stepSize;
    params.max_steps = maxSteps;
    params.max_length = 1000.0f;
    params.bounds = geomlib_trimesh_compute_aabb(&mesh.mesh);

    geomlib_streamline streamline = geomlib_streamline_create(maxSteps);

    int result = geomlib_trace_streamline(
        &mesh.mesh,
        &velocity.field,
        seed,
        &params,
        &streamline
    );

    val output = val::object();

    if (result == GEOMLIB_MESH_SUCCESS) {
        output.set("success", true);
        output.set("pointCount", (int)streamline.point_count);

        val points = val::array();
        for (size_t i = 0; i < streamline.point_count; i++) {
            points.call<void>("push", streamline.points[i].x);
            points.call<void>("push", streamline.points[i].y);
            points.call<void>("push", streamline.points[i].z);
        }
        output.set("points", points);
    } else {
        output.set("success", false);
    }

    geomlib_streamline_free(&streamline);
    return output;
}

// ============================================================================
// Tetrahedral Mesh (CFD)
// ============================================================================

class TetraMesh {
public:
    geomlib_tetramesh mesh;

    TetraMesh() {
        geomlib_tetramesh_init(&mesh, 64, 64);
    }

    ~TetraMesh() {
        geomlib_tetramesh_free(&mesh);
    }

    TetraMesh(const TetraMesh&) = delete;
    TetraMesh& operator=(const TetraMesh&) = delete;

    int getVertexCount() const {
        return mesh.vertex_count;
    }

    int getTetrahedronCount() const {
        return mesh.tetrahedron_count;
    }

    val getVertices() const {
        val vertices = val::array();
        for (size_t i = 0; i < mesh.vertex_count; i++) {
            vertices.call<void>("push", mesh.vertices[i].x);
            vertices.call<void>("push", mesh.vertices[i].y);
            vertices.call<void>("push", mesh.vertices[i].z);
        }
        return vertices;
    }

    val getTetrahedra() const {
        val tets = val::array();
        for (size_t i = 0; i < mesh.tetrahedron_count; i++) {
            tets.call<void>("push", mesh.tetrahedra[i].vertices[0]);
            tets.call<void>("push", mesh.tetrahedra[i].vertices[1]);
            tets.call<void>("push", mesh.tetrahedra[i].vertices[2]);
            tets.call<void>("push", mesh.tetrahedra[i].vertices[3]);
        }
        return tets;
    }
};

val computeCFDQuality(const TetraMesh& tetramesh) {
    geomlib_cfd_quality_thresholds thresholds = geomlib_cfd_quality_thresholds_default();

    size_t failed_count = 0;
    bool all_valid = geomlib_check_cfd_mesh_quality(
        &tetramesh.mesh,
        &thresholds,
        &failed_count
    );

    val output = val::object();
    output.set("allValid", all_valid);
    output.set("failedCount", (unsigned int)failed_count);
    output.set("passRate", (float)(tetramesh.mesh.tetrahedron_count - failed_count) /
        (float)tetramesh.mesh.tetrahedron_count * 100.0f);

    return output;
}

TetraMesh* tetrahedralizePoints(val points_array) {
    // Convert JavaScript array to C array
    size_t count = points_array["length"].as<size_t>() / 3;
    geomlib_vec3 *points = (geomlib_vec3*)malloc(count * sizeof(geomlib_vec3));

    for (size_t i = 0; i < count; i++) {
        points[i].x = points_array[i * 3 + 0].as<float>();
        points[i].y = points_array[i * 3 + 1].as<float>();
        points[i].z = points_array[i * 3 + 2].as<float>();
    }

    TetraMesh* tetramesh = new TetraMesh();
    int result = geomlib_tetrahedralize_points(points, count, &tetramesh->mesh);

    free(points);

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete tetramesh;
        return nullptr;
    }

    return tetramesh;
}

TetraMesh* tetrahedralizeBox(
    float minX, float minY, float minZ,
    float maxX, float maxY, float maxZ,
    float elementSize
) {
    geomlib_vec3 min_corner = geomlib_vec3_create(minX, minY, minZ);
    geomlib_vec3 max_corner = geomlib_vec3_create(maxX, maxY, maxZ);

    geomlib_tetra_params params = geomlib_tetra_params_default(GEOMLIB_TETRA_DELAUNAY);
    params.element_size = elementSize;

    TetraMesh* tetramesh = new TetraMesh();
    int result = geomlib_tetrahedralize_box(min_corner, max_corner, &params, &tetramesh->mesh);

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete tetramesh;
        return nullptr;
    }

    return tetramesh;
}

TetraMesh* tetrahedralizeFromSurface(const Mesh& surface, float elementSize) {
    geomlib_tetra_params params = geomlib_tetra_params_default(GEOMLIB_TETRA_DELAUNAY);
    params.element_size = elementSize;

    TetraMesh* tetramesh = new TetraMesh();
    int result = geomlib_tetrahedralize_surface(&surface.mesh, &params, &tetramesh->mesh);

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete tetramesh;
        return nullptr;
    }

    return tetramesh;
}

val computeTetraMeshStats(const TetraMesh& tetramesh) {
    geomlib_tetramesh_stats stats;
    geomlib_tetramesh_compute_stats(&tetramesh.mesh, &stats);

    val output = val::object();
    output.set("vertexCount", (unsigned int)stats.vertex_count);
    output.set("tetrahedronCount", (unsigned int)stats.tetrahedron_count);
    output.set("boundaryFaceCount", (unsigned int)stats.boundary_face_count);
    output.set("totalVolume", stats.total_volume);
    output.set("minVolume", stats.min_volume);
    output.set("maxVolume", stats.max_volume);
    output.set("avgVolume", stats.avg_volume);
    output.set("minQuality", stats.min_quality);
    output.set("avgQuality", stats.avg_quality);
    output.set("worstAspectRatio", stats.worst_aspect_ratio);
    output.set("sliverCount", (unsigned int)stats.sliver_count);
    output.set("boundaryVertices", (unsigned int)stats.boundary_vertices);

    return output;
}

Mesh* extractBoundaryFromTetraMesh(const TetraMesh& tetramesh) {
    Mesh* boundary = new Mesh();
    int result = geomlib_tetramesh_extract_boundary(&tetramesh.mesh, &boundary->mesh);

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete boundary;
        return nullptr;
    }

    return boundary;
}

// ============================================================================
// Subdivision Surfaces
// ============================================================================

Mesh* subdivideLoop(const Mesh& input, int levels) {
    Mesh* output = new Mesh();

    geomlib_subdivision_params params = geomlib_subdivision_params_default(
        GEOMLIB_SUBDIV_LOOP, levels);

    int result = geomlib_subdivide_mesh(&input.mesh, &params, &output->mesh);

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete output;
        return nullptr;
    }

    return output;
}

Mesh* subdivideMidpoint(const Mesh& input, int levels) {
    Mesh* output = new Mesh();

    geomlib_subdivision_params params = geomlib_subdivision_params_default(
        GEOMLIB_SUBDIV_MIDPOINT, levels);

    int result = geomlib_subdivide_mesh(&input.mesh, &params, &output->mesh);

    if (result != GEOMLIB_MESH_SUCCESS) {
        delete output;
        return nullptr;
    }

    return output;
}

// ============================================================================
// Mesh Optimization
// ============================================================================

Mesh* optimizeVertexCache(Mesh* mesh) {
    if (!mesh) return nullptr;

    geomlib_optimize_vertex_cache(&mesh->mesh, GEOMLIB_VCACHE_FORSYTH, 16);
    return mesh;
}

Mesh* weldVertices(Mesh* mesh, float tolerance) {
    if (!mesh) return nullptr;

    geomlib_weld_result weld = geomlib_weld_vertices(&mesh->mesh, tolerance);
    geomlib_apply_weld_result(&mesh->mesh, &weld);
    geomlib_weld_result_free(&weld);

    return mesh;
}

Mesh* optimizeVertexFetch(Mesh* mesh) {
    if (!mesh) return nullptr;

    geomlib_optimize_vertex_fetch(&mesh->mesh);
    return mesh;
}

val optimizeMeshComplete(Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_optimize_params params = geomlib_optimize_params_default();
    geomlib_optimize_result result = geomlib_optimize_mesh(&mesh->mesh, &params);

    val output = val::object();
    output.set("verticesBefore", (int)result.vertices_before);
    output.set("verticesAfter", (int)result.vertices_after);
    output.set("acmrBefore", result.cache_before.acmr);
    output.set("acmrAfter", result.cache_after.acmr);
    output.set("improvement", (result.cache_before.acmr - result.cache_after.acmr) / result.cache_before.acmr * 100.0f);
    output.set("optimizationTimeMs", result.optimization_time_ms);

    return output;
}

val analyzeCacheEfficiency(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_cache_stats stats = geomlib_analyze_vertex_cache(&mesh->mesh, 16);

    val output = val::object();
    output.set("acmr", stats.acmr);
    output.set("atvr", stats.atvr);
    output.set("cacheMisses", (int)stats.cache_misses);
    output.set("transformedVertices", (int)stats.transformed_vertices);
    output.set("cacheSize", (int)stats.cache_size);

    return output;
}

// ============================================================================
// Topology Analysis
// ============================================================================

val analyzeTopology(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_euler_characteristic euler = geomlib_compute_euler_characteristic(&mesh->mesh);

    val output = val::object();
    output.set("vertices", euler.vertices);
    output.set("edges", euler.edges);
    output.set("faces", euler.faces);
    output.set("chi", euler.chi);
    output.set("genus", euler.genus);
    output.set("boundaries", euler.boundaries);
    output.set("components", euler.components);

    return output;
}

val findConnectedComponents(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_component_info info = geomlib_find_connected_components(&mesh->mesh);

    val output = val::object();
    output.set("componentCount", (int)info.component_count);

    val vertexLabels = val::array();
    for (size_t i = 0; i < mesh->mesh.vertex_count; i++) {
        vertexLabels.call<void>("push", (int)info.vertex_labels[i]);
    }
    output.set("vertexLabels", vertexLabels);

    val componentSizes = val::array();
    for (uint32_t i = 0; i < info.component_count; i++) {
        val component = val::object();
        component.set("vertices", (int)info.component_vertex_counts[i]);
        component.set("faces", (int)info.component_face_counts[i]);
        componentSizes.call<void>("push", component);
    }
    output.set("components", componentSizes);

    geomlib_component_info_free(&info);

    return output;
}

val extractBoundaryLoops(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_boundary_loops loops = geomlib_extract_boundary_loops(&mesh->mesh);

    val output = val::object();
    output.set("loopCount", (int)loops.loop_count);

    val loopArray = val::array();
    for (size_t i = 0; i < loops.loop_count; i++) {
        val loop = val::object();
        loop.set("vertexCount", (int)loops.loops[i].vertex_count);
        loop.set("length", loops.loops[i].length);

        val vertices = val::array();
        for (size_t j = 0; j < loops.loops[i].vertex_count; j++) {
            vertices.call<void>("push", (int)loops.loops[i].vertices[j]);
        }
        loop.set("vertices", vertices);

        loopArray.call<void>("push", loop);
    }
    output.set("loops", loopArray);

    geomlib_boundary_loops_free(&loops);

    return output;
}

val analyzeManifold(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&mesh->mesh);

    val output = val::object();
    output.set("isManifold", analysis.is_manifold);
    output.set("isEdgeManifold", analysis.is_edge_manifold);
    output.set("isVertexManifold", analysis.is_vertex_manifold);
    output.set("isOriented", analysis.is_oriented);
    output.set("isClosed", analysis.is_closed);
    output.set("nonmanifoldCount", (int)analysis.nonmanifold_count);

    geomlib_manifold_analysis_free(&analysis);

    return output;
}

val computeValenceStats(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_valence_stats stats = geomlib_compute_valence_stats(&mesh->mesh);

    val output = val::object();
    output.set("minValence", (int)stats.min_valence);
    output.set("maxValence", (int)stats.max_valence);
    output.set("avgValence", stats.avg_valence);

    val histogram = val::array();
    for (size_t i = 0; i < stats.histogram_size; i++) {
        histogram.call<void>("push", (int)stats.histogram[i]);
    }
    output.set("histogram", histogram);

    geomlib_valence_stats_free(&stats);

    return output;
}

// ============================================================================
// Hole Filling
// ============================================================================

val detectHoles(const Mesh* mesh) {
    if (!mesh) {
        return val::object();
    }

    geomlib_holes holes = geomlib_detect_holes(&mesh->mesh);

    val output = val::object();
    output.set("holeCount", (int)holes.hole_count);

    val holesArray = val::array();
    for (size_t i = 0; i < holes.hole_count; i++) {
        val hole = val::object();
        hole.set("vertexCount", (int)holes.holes[i].vertex_count);
        hole.set("perimeter", holes.holes[i].perimeter);
        hole.set("areaEstimate", holes.holes[i].area_estimate);
        hole.set("isPlanar", holes.holes[i].is_planar);
        hole.set("planarity", holes.holes[i].planarity);

        val centroid = val::object();
        centroid.set("x", holes.holes[i].centroid.x);
        centroid.set("y", holes.holes[i].centroid.y);
        centroid.set("z", holes.holes[i].centroid.z);
        hole.set("centroid", centroid);

        val normal = val::object();
        normal.set("x", holes.holes[i].normal.x);
        normal.set("y", holes.holes[i].normal.y);
        normal.set("z", holes.holes[i].normal.z);
        hole.set("normal", normal);

        // Boundary vertices
        val boundaryVertices = val::array();
        for (size_t j = 0; j < holes.holes[i].vertex_count; j++) {
            boundaryVertices.set(j, holes.holes[i].boundary_vertices[j]);
        }
        hole.set("boundaryVertices", boundaryVertices);

        holesArray.set(i, hole);
    }
    output.set("holes", holesArray);

    geomlib_holes_free(&holes);

    return output;
}

Mesh* fillAllHoles(Mesh* mesh, const std::string& method) {
    if (!mesh) {
        return nullptr;
    }

    geomlib_fill_method fill_method = GEOMLIB_FILL_TANGENT_PLANE;
    if (method == "advancing_front") {
        fill_method = GEOMLIB_FILL_ADVANCING_FRONT;
    } else if (method == "minimal_surface") {
        fill_method = GEOMLIB_FILL_MINIMAL_SURFACE;
    } else if (method == "natural") {
        fill_method = GEOMLIB_FILL_NATURAL;
    }

    geomlib_fill_params params = geomlib_fill_params_default(fill_method);

    int result = geomlib_fill_all_holes(&mesh->mesh, &params);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return nullptr;
    }

    return mesh;
}

// ============================================================================
// Feature Detection
// ============================================================================

val detectFeatures(const Mesh* mesh, float sharpAngleThreshold) {
    if (!mesh) {
        return val::object();
    }

    geomlib_feature_params params = geomlib_feature_params_default();
    params.sharp_angle_threshold = sharpAngleThreshold;

    geomlib_features features = geomlib_detect_features(&mesh->mesh, &params);

    val output = val::object();
    output.set("edgeCount", (int)features.edge_count);
    output.set("vertexCount", (int)features.vertex_count);
    output.set("curveCount", (int)features.curve_count);

    // Feature edges
    val edges = val::array();
    for (size_t i = 0; i < features.edge_count; i++) {
        val edge = val::object();
        edge.set("v0", features.edges[i].v0);
        edge.set("v1", features.edges[i].v1);
        edge.set("type", (int)features.edges[i].type);
        edge.set("sharpness", features.edges[i].sharpness);

        std::string type_name;
        switch (features.edges[i].type) {
            case GEOMLIB_FEATURE_SHARP_EDGE: type_name = "sharp_edge"; break;
            case GEOMLIB_FEATURE_RIDGE: type_name = "ridge"; break;
            case GEOMLIB_FEATURE_VALLEY: type_name = "valley"; break;
            case GEOMLIB_FEATURE_BOUNDARY: type_name = "boundary"; break;
            case GEOMLIB_FEATURE_CORNER: type_name = "corner"; break;
            default: type_name = "none"; break;
        }
        edge.set("typeName", type_name);

        edges.set(i, edge);
    }
    output.set("edges", edges);

    // Feature vertices (corners)
    val vertices = val::array();
    for (size_t i = 0; i < features.vertex_count; i++) {
        val vertex = val::object();
        vertex.set("index", features.vertices[i].vertex_index);
        vertex.set("edgeCount", features.vertices[i].edge_count);
        vertex.set("sharpness", features.vertices[i].sharpness);
        vertices.set(i, vertex);
    }
    output.set("vertices", vertices);

    // Statistics
    geomlib_feature_stats stats = geomlib_compute_feature_stats(&features);
    val statsObj = val::object();
    statsObj.set("sharpEdgeCount", (int)stats.sharp_edge_count);
    statsObj.set("ridgeCount", (int)stats.ridge_count);
    statsObj.set("valleyCount", (int)stats.valley_count);
    statsObj.set("cornerCount", (int)stats.corner_count);
    statsObj.set("boundaryEdgeCount", (int)stats.boundary_edge_count);
    statsObj.set("avgSharpness", stats.avg_sharpness);
    statsObj.set("maxSharpness", stats.max_sharpness);
    statsObj.set("minSharpness", stats.min_sharpness);
    output.set("stats", statsObj);

    geomlib_features_free(&features);

    return output;
}

val detectSharpEdges(const Mesh* mesh, float angleThreshold) {
    if (!mesh) {
        return val::object();
    }

    size_t edge_count;
    geomlib_feature_edge *edges = geomlib_detect_sharp_edges(
        &mesh->mesh, angleThreshold, &edge_count);

    val output = val::object();
    output.set("edgeCount", (int)edge_count);

    val edgesArray = val::array();
    for (size_t i = 0; i < edge_count; i++) {
        val edge = val::object();
        edge.set("v0", edges[i].v0);
        edge.set("v1", edges[i].v1);
        edge.set("sharpness", edges[i].sharpness);

        edgesArray.set(i, edge);
    }
    output.set("edges", edgesArray);

    free(edges);

    return output;
}

// ============================================================================
// SDF (Signed Distance Field)
// ============================================================================

class SDF {
public:
    geomlib_sdf sdf;

    SDF() {
        memset(&sdf, 0, sizeof(geomlib_sdf));
    }

    ~SDF() {
        geomlib_sdf_free(&sdf);
    }

    SDF(const SDF&) = delete;
    SDF& operator=(const SDF&) = delete;

    int getResolutionX() const { return sdf.resolution[0]; }
    int getResolutionY() const { return sdf.resolution[1]; }
    int getResolutionZ() const { return sdf.resolution[2]; }
    float getVoxelSize() const { return sdf.voxel_size; }

    float queryDistance(float x, float y, float z) const {
        geomlib_vec3 point = geomlib_vec3_create(x, y, z);
        return geomlib_sdf_query_distance(&sdf, point);
    }

    bool isInside(float x, float y, float z) const {
        geomlib_vec3 point = geomlib_vec3_create(x, y, z);
        return geomlib_sdf_is_inside(&sdf, point);
    }

    val getStats() const {
        geomlib_sdf_stats stats = geomlib_compute_sdf_stats(&sdf);
        val output = val::object();
        output.set("minDistance", stats.min_distance);
        output.set("maxDistance", stats.max_distance);
        output.set("avgDistance", stats.avg_distance);
        output.set("insideVoxels", (int)stats.inside_voxels);
        output.set("outsideVoxels", (int)stats.outside_voxels);
        output.set("surfaceVoxels", (int)stats.surface_voxels);
        output.set("memoryBytes", (int)stats.memory_bytes);
        return output;
    }
};

SDF* generateSDF(const Mesh* mesh, float voxelSize) {
    if (!mesh) {
        return nullptr;
    }

    SDF* sdf = new SDF();
    sdf->sdf = geomlib_mesh_to_sdf_voxel_size(&mesh->mesh, voxelSize);

    return sdf;
}

SDF* generateSDFResolution(const Mesh* mesh, int resolution) {
    if (!mesh) {
        return nullptr;
    }

    SDF* sdf = new SDF();
    sdf->sdf = geomlib_mesh_to_sdf_resolution(&mesh->mesh, resolution);

    return sdf;
}

// ============================================================================
// Mesh Statistics
// ============================================================================

val computeMeshStatistics(const Mesh* mesh, int histogramBins) {
    if (!mesh) {
        return val::object();
    }

    geomlib_mesh_stats stats = geomlib_compute_mesh_stats(&mesh->mesh, histogramBins);

    val output = val::object();
    output.set("vertexCount", (int)stats.vertex_count);
    output.set("triangleCount", (int)stats.triangle_count);
    output.set("edgeCount", (int)stats.edge_count);
    output.set("isClosed", stats.is_closed);
    output.set("surfaceArea", stats.surface_area);
    output.set("volume", stats.volume);

    // Edge stats
    val edges = val::object();
    edges.set("min", stats.edges.min_length);
    edges.set("max", stats.edges.max_length);
    edges.set("avg", stats.edges.avg_length);
    edges.set("median", stats.edges.median_length);
    edges.set("stdDev", stats.edges.std_dev);
    output.set("edges", edges);

    // Angle stats
    val angles = val::object();
    angles.set("min", stats.angles.min_angle_deg);
    angles.set("max", stats.angles.max_angle_deg);
    angles.set("avg", stats.angles.avg_angle_deg);
    angles.set("degenerateCount", (int)stats.angles.degenerate_count);
    angles.set("poorQualityCount", (int)stats.angles.poor_quality_count);
    output.set("angles", angles);

    // Area stats
    val areas = val::object();
    areas.set("min", stats.areas.min_area);
    areas.set("max", stats.areas.max_area);
    areas.set("avg", stats.areas.avg_area);
    areas.set("total", stats.areas.total_area);
    areas.set("stdDev", stats.areas.std_dev);
    output.set("areas", areas);

    // Aspect ratio stats
    val aspectRatios = val::object();
    aspectRatios.set("min", stats.aspect_ratios.min_ratio);
    aspectRatios.set("max", stats.aspect_ratios.max_ratio);
    aspectRatios.set("avg", stats.aspect_ratios.avg_ratio);
    aspectRatios.set("goodCount", (int)stats.aspect_ratios.good_count);
    aspectRatios.set("acceptableCount", (int)stats.aspect_ratios.acceptable_count);
    aspectRatios.set("poorCount", (int)stats.aspect_ratios.poor_count);
    output.set("aspectRatios", aspectRatios);

    geomlib_mesh_stats_free(&stats);

    return output;
}

// ============================================================================
// Embind Bindings
// ============================================================================

EMSCRIPTEN_BINDINGS(geomlib) {
    // Vec3
    class_<Vec3>("Vec3")
        .constructor<>()
        .constructor<float, float, float>()
        .property("x", &Vec3::x)
        .property("y", &Vec3::y)
        .property("z", &Vec3::z)
        ;

    // Mesh
    class_<Mesh>("Mesh")
        .constructor<>()
        .function("addVertex", &Mesh::addVertex)
        .function("addTriangle", &Mesh::addTriangle)
        .function("getVertexCount", &Mesh::getVertexCount)
        .function("getTriangleCount", &Mesh::getTriangleCount)
        .function("getVertices", &Mesh::getVertices)
        .function("getIndices", &Mesh::getIndices)
        .function("getNormals", &Mesh::getNormals)
        .function("getSurfaceArea", &Mesh::getSurfaceArea)
        .function("getVolume", &Mesh::getVolume)
        .function("getCentroid", &Mesh::getCentroid)
        ;

    // BVH
    class_<BVH>("BVH")
        .constructor<const Mesh*>()
        .function("traceRay", &BVH::traceRay)
        .function("getStats", &BVH::getStats)
        .function("getBoundingBoxes", &BVH::getBoundingBoxes, allow_raw_pointers())
        ;

    // Mesh creation functions
    function("createBox", &createBox, allow_raw_pointers());
    function("createSphere", &createSphere, allow_raw_pointers());
    function("createCylinder", &createCylinder, allow_raw_pointers());
    function("createPlane", &createPlane, allow_raw_pointers());

    // Collision proxy generation
    function("generateCollisionProxy", &generateCollisionProxy, allow_raw_pointers());
    function("simplifyMesh", &simplifyMesh, allow_raw_pointers());
    function("computeConvexHull", &computeConvexHull, allow_raw_pointers());

    // BVH
    function("buildBVH", &buildBVH, allow_raw_pointers());

    // Analysis
    function("analyzeMeshQuality", &analyzeMeshQuality);
    function("estimateCollisionCost", &estimateCollisionCost);

    // Point Cloud
    class_<PointCloud>("PointCloud")
        .constructor<>()
        .function("addPoint", &PointCloud::addPoint)
        .function("addPointWithNormal", &PointCloud::addPointWithNormal)
        .function("getPointCount", &PointCloud::getPointCount)
        .function("getPoints", &PointCloud::getPoints)
        .function("getStats", &PointCloud::getStats)
        ;

    // Surface reconstruction
    function("gaussiansToPointCloud", &gaussiansToPointCloud, allow_raw_pointers());
    function("gaussiansSamplePoints", &gaussiansSamplePoints, allow_raw_pointers());
    function("reconstructSurface", &reconstructSurface, allow_raw_pointers());
    function("reconstructSurfaceAuto", &reconstructSurfaceAuto, allow_raw_pointers());

    // Mesh healing
    function("fixFaceOrientation", &fixFaceOrientation, allow_raw_pointers());
    function("removeDegenerateTriangles", &removeDegenerateTriangles, allow_raw_pointers());
    function("removeDuplicateVertices", &removeDuplicateVertices, allow_raw_pointers());
    function("repairMesh", &repairMesh, allow_raw_pointers());
    function("analyzeMeshTopology", &analyzeMeshTopology);

    // Remeshing
    function("remeshUniform", &remeshUniform, allow_raw_pointers());
    function("remeshAdaptive", &remeshAdaptive, allow_raw_pointers());
    function("computeTargetEdgeLength", &computeTargetEdgeLength);
    function("computeMeanEdgeLength", &computeMeanEdgeLength);
    function("getRemeshingStats", &getRemeshingStats);

    // Parameterization
    function("parameterizeHarmonic", &parameterizeHarmonic);
    function("parameterizeLSCM", &parameterizeLSCM);

    // Curvature
    function("computeCurvature", &computeCurvature);
    function("computeMeanCurvature", &computeMeanCurvature);
    function("computeGaussianCurvature", &computeGaussianCurvature);

    // Geodesic Distances
    function("computeGeodesicDistances", &computeGeodesicDistances);
    function("computeGeodesicDistancesMultiple", &computeGeodesicDistancesMultiple);

    // Decimation
    function("decimateQEM", &decimateQEM);
    function("decimateToVertexCount", &decimateToVertexCount);

    // Smoothing
    function("smoothLaplacian", &smoothLaplacian);
    function("smoothTaubin", &smoothTaubin);
    function("smoothNormals", &smoothNormals);

    // Mesh Quality
    function("computeMeshQuality", &computeMeshQuality);
    function("checkMeshQuality", &checkMeshQuality);

    // Texture Atlas
    function("generateTextureAtlas", &generateTextureAtlas);
    function("generateTextureAtlasCustom", &generateTextureAtlasCustom);

    // CFD Utilities
    class_<ScalarField>("ScalarField")
        .constructor<int>()
        .function("setValue", &ScalarField::setValue)
        .function("getValue", &ScalarField::getValue)
        .function("getValues", &ScalarField::getValues)
        ;

    class_<VectorField>("VectorField")
        .constructor<int>()
        .function("setVector", &VectorField::setVector)
        .function("getVectors", &VectorField::getVectors)
        ;

    function("computeGradient", &computeGradient, allow_raw_pointers());
    function("traceStreamline", &traceStreamline);
    function("computeCFDQuality", &computeCFDQuality);

    // Tetrahedral Mesh (CFD)
    class_<TetraMesh>("TetraMesh")
        .constructor<>()
        .function("getVertexCount", &TetraMesh::getVertexCount)
        .function("getTetrahedronCount", &TetraMesh::getTetrahedronCount)
        .function("getVertices", &TetraMesh::getVertices)
        .function("getTetrahedra", &TetraMesh::getTetrahedra)
        ;

    function("tetrahedralizePoints", &tetrahedralizePoints, allow_raw_pointers());
    function("tetrahedralizeBox", &tetrahedralizeBox, allow_raw_pointers());
    function("tetrahedralizeFromSurface", &tetrahedralizeFromSurface, allow_raw_pointers());
    function("computeTetraMeshStats", &computeTetraMeshStats);
    function("extractBoundaryFromTetraMesh", &extractBoundaryFromTetraMesh, allow_raw_pointers());

    // Mesh Optimization
    function("optimizeVertexCache", select_overload<Mesh*(Mesh*)>(&optimizeVertexCache), allow_raw_pointers());
    function("weldVertices", select_overload<Mesh*(Mesh*, float)>(&weldVertices), allow_raw_pointers());
    function("optimizeVertexFetch", select_overload<Mesh*(Mesh*)>(&optimizeVertexFetch), allow_raw_pointers());
    function("optimizeMeshComplete", select_overload<val(Mesh*)>(&optimizeMeshComplete), allow_raw_pointers());
    function("analyzeCacheEfficiency", select_overload<val(const Mesh*)>(&analyzeCacheEfficiency), allow_raw_pointers());

    // Topology Analysis
    function("analyzeTopology", select_overload<val(const Mesh*)>(&analyzeTopology), allow_raw_pointers());
    function("findConnectedComponents", select_overload<val(const Mesh*)>(&findConnectedComponents), allow_raw_pointers());
    function("extractBoundaryLoops", select_overload<val(const Mesh*)>(&extractBoundaryLoops), allow_raw_pointers());
    function("analyzeManifold", select_overload<val(const Mesh*)>(&analyzeManifold), allow_raw_pointers());
    function("computeValenceStats", select_overload<val(const Mesh*)>(&computeValenceStats), allow_raw_pointers());

    // Mesh Statistics
    function("computeMeshStatistics", select_overload<val(const Mesh*, int)>(&computeMeshStatistics), allow_raw_pointers());

    // Hole Filling
    function("detectHoles", select_overload<val(const Mesh*)>(&detectHoles), allow_raw_pointers());
    function("fillAllHoles", select_overload<Mesh*(Mesh*, const std::string&)>(&fillAllHoles), allow_raw_pointers());

    // Feature Detection
    function("detectFeatures", select_overload<val(const Mesh*, float)>(&detectFeatures), allow_raw_pointers());
    function("detectSharpEdges", select_overload<val(const Mesh*, float)>(&detectSharpEdges), allow_raw_pointers());

    // SDF (Signed Distance Field)
    class_<SDF>("SDF")
        .constructor<>()
        .function("getResolutionX", &SDF::getResolutionX)
        .function("getResolutionY", &SDF::getResolutionY)
        .function("getResolutionZ", &SDF::getResolutionZ)
        .function("getVoxelSize", &SDF::getVoxelSize)
        .function("queryDistance", &SDF::queryDistance)
        .function("isInside", &SDF::isInside)
        .function("getStats", &SDF::getStats)
        ;

    function("generateSDF", &generateSDF, allow_raw_pointers());
    function("generateSDFResolution", &generateSDFResolution, allow_raw_pointers());

    // Subdivision
    function("subdivideLoop", &subdivideLoop, allow_raw_pointers());
    function("subdivideMidpoint", &subdivideMidpoint, allow_raw_pointers());
}
