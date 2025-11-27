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

    // Subdivision
    function("subdivideLoop", &subdivideLoop, allow_raw_pointers());
    function("subdivideMidpoint", &subdivideMidpoint, allow_raw_pointers());
}
