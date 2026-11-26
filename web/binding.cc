/*
 * geomlib - WebAssembly/JavaScript Bindings
 * Using Emscripten embind
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <string>
#include <cmath>
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

    val getUVs() const {
        val uvs = val::array();
        if (mesh.texcoords) {
            for (size_t i = 0; i < mesh.vertex_count; i++) {
                uvs.call<void>("push", mesh.texcoords[i].x);
                uvs.call<void>("push", mesh.texcoords[i].y);
            }
        }
        return uvs;
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

    geomlib_surface_params params = geomlib_surface_params_default(
        (geomlib_surface_method)method);
    params.voxel_size = voxelSize;

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
    VectorField* gradient = new VectorField(0);
    geomlib_vector_field_free(&gradient->field);

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
    size_t cells = tetramesh.mesh.tetrahedron_count;
    output.set("passRate", cells ? 100.0f * (cells - failed_count) / cells : 100.0f);

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

Mesh* subdivideButterfly(const Mesh& input, int levels) {
    Mesh* output = new Mesh();

    geomlib_subdivision_params params = geomlib_subdivision_params_default(
        GEOMLIB_SUBDIV_BUTTERFLY, levels);

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

    Mesh* extractSurface(float isoValue) const {
        geomlib_trimesh extracted = {};
        if (geomlib_sdf_to_mesh(&sdf, isoValue, &extracted) != GEOMLIB_MESH_SUCCESS)
            return nullptr;

        Mesh* result = new Mesh();
        geomlib_trimesh_free(&result->mesh);
        result->mesh = extracted;
        return result;
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
    if (!sdf->sdf.values) {
        delete sdf;
        return nullptr;
    }
    return sdf;
}

SDF* generateSDFResolution(const Mesh* mesh, int resolution) {
    if (!mesh || resolution <= 0) {
        return nullptr;
    }

    SDF* sdf = new SDF();
    sdf->sdf = geomlib_mesh_to_sdf_resolution(&mesh->mesh, resolution);
    if (!sdf->sdf.values) {
        delete sdf;
        return nullptr;
    }
    return sdf;
}

// ============================================================================
// Mesh Statistics
// ============================================================================

val computeMeshStatistics(const Mesh* mesh, int histogramBins) {
    if (!mesh) {
        return val::object();
    }

    geomlib_mesh_stats stats = geomlib_compute_mesh_stats(&mesh->mesh, histogramBins > 0 ? histogramBins : 0);

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

// ============================================================================
// Wave 1: parameterization and Catmull-Clark
// ============================================================================

bool parameterizeMesh(Mesh* mesh, int method, int boundary, int iterations) {
    if (!mesh || !mesh->mesh.vertices) return false;
    geomlib_uv_coords uv = geomlib_uv_coords_create(mesh->mesh.vertex_count,
                                                    mesh->mesh.triangle_count);
    geomlib_param_params params = geomlib_param_params_default(
        static_cast<geomlib_param_method>(method));
    params.boundary = static_cast<geomlib_boundary_shape>(boundary);
    params.max_iterations = iterations;
    int rc = geomlib_parameterize(&mesh->mesh, &params, &uv);
    if (rc == GEOMLIB_MESH_SUCCESS) {
        free(mesh->mesh.texcoords);
        mesh->mesh.texcoords = (geomlib_vec2*)malloc(
            mesh->mesh.vertex_count * sizeof(geomlib_vec2));
        if (mesh->mesh.texcoords)
            memcpy(mesh->mesh.texcoords, uv.uvs,
                   mesh->mesh.vertex_count * sizeof(geomlib_vec2));
    }
    geomlib_uv_coords_free(&uv);
    return rc == GEOMLIB_MESH_SUCCESS;
}

val reconstructRBF(val points) {
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud, 512);
    size_t n = points["length"].as<size_t>() / 3;
    for (size_t i = 0; i < n; ++i) {
        geomlib_vec3 p = geomlib_vec3_create(points[i*3].as<float>(),
                                             points[i*3+1].as<float>(),
                                             points[i*3+2].as<float>());
        geomlib_pointcloud_add_point(&cloud, p, nullptr, 1.0f);
    }
    geomlib_rbf_params params = geomlib_rbf_params_default();
    params.resolution = 20;
    params.max_points = (int)n;
    Mesh* out = new Mesh();
    int rc = geomlib_reconstruct_rbf(&cloud, &params, &out->mesh);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("triangles", static_cast<int>(out->mesh.triangle_count));
    delete out;
    geomlib_pointcloud_free(&cloud);
    return result;
}

val metricDemo() {
    geomlib_metric3 M=geomlib_metric_diagonal(4.0,1.0,0.25);
    double d=geomlib_metric_distance(M,geomlib_vec3_create(0,0,0),
                                     geomlib_vec3_create(2,3,6));
    double a=geomlib_metric_anisotropy(M);
    val result=val::object();
    result.set("success",std::fabs(d-std::sqrt(34.0))<1e-4 &&
                         std::fabs(a-4.0)<1e-4);
    result.set("distance",d);
    return result;
}

val geodesicMmpDemo() {
    geomlib_trimesh cube;
    geomlib_trimesh_init(&cube, 8, 12);
    geomlib_vec3 v[8]={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(1,1,0),geomlib_vec3_create(0,1,0),
        geomlib_vec3_create(0,0,1),geomlib_vec3_create(1,0,1),
        geomlib_vec3_create(1,1,1),geomlib_vec3_create(0,1,1)};
    for (int i=0;i<8;++i) geomlib_trimesh_add_vertex(&cube,v[i],nullptr,nullptr);
    int f[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},
                  {3,2,6},{3,6,7},{0,3,7},{0,7,4},{1,5,6},{1,6,2}};
    for (int i=0;i<12;++i)
        geomlib_trimesh_add_triangle(&cube,f[i][0],f[i][1],f[i][2]);
    double len=0.0;
    int rc=geomlib_geodesic_mmp(&cube,0,6,&len);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(len-2.2360679)<1e-4);
    result.set("length",len);
    geomlib_trimesh_free(&cube);
    return result;
}

val circleEllipseDemo() {
    geomlib_circle_ellipse_result r;
    int hit=geomlib_circle_ellipse_overlap(
        geomlib_vec2_create(3,0),1.5f,geomlib_vec2_create(0,0),
        2.0f,1.0f,0.0f,&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.5f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val polygonIntersectionAreaDemo() {
    std::vector<geomlib_vec2> a={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)};
    std::vector<geomlib_vec2> b={
        geomlib_vec2_create(0.5f,0.5f),geomlib_vec2_create(1.5f,0.5f),
        geomlib_vec2_create(1.5f,1.5f),geomlib_vec2_create(0.5f,1.5f)};
    double area=geomlib_polygon_intersection_area(a.data(),a.size(),
                                                  b.data(),b.size());
    val result=val::object();
    result.set("success",std::fabs(area-0.25)<1e-3);
    result.set("area",area);
    return result;
}

val circleAreaDemo() {
    double a=geomlib_circle_intersection_area(1.0,1.0,1.0);
    val result=val::object();
    result.set("success",std::fabs(a-1.228370)<1e-4);
    result.set("area",a);
    return result;
}

val arcLengthDemo() {
    std::vector<geomlib_vec3> p={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(3,0,0),
        geomlib_vec3_create(3,4,0)};
    double total=geomlib_polyline_length(p.data(),p.size());
    val result=val::object();
    result.set("success",std::fabs(total-7.0)<1e-4);
    result.set("length",total);
    return result;
}

val line2dDemo() {
    geomlib_line2d_intersection r=geomlib_line2d_intersect(
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(1,-1),geomlib_vec2_create(1,1));
    val result=val::object();
    result.set("success",r.kind==GEOMLIB_LINE2D_POINT &&
                         std::fabs(r.point.x-1.0f)<1e-3 &&
                         std::fabs(r.point.y)<1e-3);
    return result;
}

val circleTriangleDemo() {
    geomlib_circle_triangle_result r;
    int hit=geomlib_circle_triangle_intersect(
        geomlib_vec2_create(0.5f,0.5f),0.1f,
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(0,2),&r);
    val result=val::object();
    result.set("success",hit && r.center_inside==1);
    return result;
}

val ellipseMeasuresDemo() {
    double p=geomlib_ellipse_perimeter(2.0f,1.0f);
    val result=val::object();
    result.set("success",std::fabs(p-9.688448)<1e-3);
    result.set("perimeter",p);
    return result;
}

val ellipseClosestDemo() {
    geomlib_ellipse_closest r;
    int ok=geomlib_ellipse_closest_point(
        geomlib_vec2_create(3,0),geomlib_vec2_create(0,0),
        2.0f,1.0f,0.0f,&r);
    val result=val::object();
    result.set("success",ok && std::fabs(r.distance-1.0f)<1e-3 &&
                         std::fabs(r.point.x-2.0f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val capsuleTriangleDemo() {
    geomlib_capsule_triangle_result r;
    int hit=geomlib_capsule_triangle_overlap(
        geomlib_vec3_create(0.5f,0.5f,-2),geomlib_vec3_create(0.5f,0.5f,2),
        0.2f,geomlib_vec3_create(0,0,0),geomlib_vec3_create(2,0,0),
        geomlib_vec3_create(0,2,0),&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.2f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val ellipseDemo() {
    double v=geomlib_ellipse_implicit(
        geomlib_vec2_create(1,0),geomlib_vec2_create(0,0),2.0f,1.0f,0.0f);
    val result=val::object();
    result.set("success",std::fabs(v-0.25)<1e-4);
    result.set("value",v);
    return result;
}

val obbObbDemo() {
    geomlib_obb a,b;
    a.center=geomlib_vec3_create(0,0,0);
    a.axes[0]=geomlib_vec3_create(1,0,0);
    a.axes[1]=geomlib_vec3_create(0,1,0);
    a.axes[2]=geomlib_vec3_create(0,0,1);
    a.half_extents=geomlib_vec3_create(1,1,1);
    b=a;
    b.center=geomlib_vec3_create(5,0,0);
    val result=val::object();
    result.set("success",geomlib_obb_obb_overlap(&a,&b)==0);
    return result;
}

val aabbObbDemo() {
    geomlib_obb box;
    box.center=geomlib_vec3_create(6,1,1);
    box.axes[0]=geomlib_vec3_create(1,0,0);
    box.axes[1]=geomlib_vec3_create(0,1,0);
    box.axes[2]=geomlib_vec3_create(0,0,1);
    box.half_extents=geomlib_vec3_create(1,1,1);
    int hit=geomlib_aabb_obb_overlap(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(2,2,2),&box);
    val result=val::object();
    result.set("success",hit==0);
    return result;
}

val triObbDemo() {
    geomlib_obb box;
    box.center=geomlib_vec3_create(0,0,0);
    box.axes[0]=geomlib_vec3_create(1,0,0);
    box.axes[1]=geomlib_vec3_create(0,1,0);
    box.axes[2]=geomlib_vec3_create(0,0,1);
    box.half_extents=geomlib_vec3_create(1,1,1);
    int hit=geomlib_triangle_obb_overlap(
        geomlib_vec3_create(0.5f,0.5f,-2),
        geomlib_vec3_create(0.5f,1.5f,0.5f),
        geomlib_vec3_create(1.5f,0.5f,0.5f),&box);
    val result=val::object();
    result.set("success",hit==1);
    return result;
}

val rayConeDemo() {
    geomlib_ray_cone_hit h;
    int hit=geomlib_ray_cone_intersect(
        geomlib_vec3_create(-3,0,2),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),
        (float)(std::acos(std::sqrt(2.0)/2.0)),&h);
    val result=val::object();
    result.set("success",hit && std::fabs(h.t-1.0)<1e-3 &&
                         std::fabs(h.point.x+2.0f)<1e-3);
    result.set("t",h.t);
    return result;
}

val capsuleObbDemo() {
    geomlib_obb box;
    box.center=geomlib_vec3_create(0,0,0);
    box.axes[0]=geomlib_vec3_create(1,0,0);
    box.axes[1]=geomlib_vec3_create(0,1,0);
    box.axes[2]=geomlib_vec3_create(0,0,1);
    box.half_extents=geomlib_vec3_create(1,1,1);
    geomlib_capsule_obb_result r;
    int hit=geomlib_capsule_obb_overlap(
        geomlib_vec3_create(2,0.5f,0.5f),geomlib_vec3_create(3,0.5f,0.5f),
        1.5f,&box,&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.5f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val circleAabbDemo() {
    geomlib_circle_aabb_result r;
    int hit=geomlib_circle_aabb_overlap(
        geomlib_vec2_create(2,0.5f),1.5f,
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,1),&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.5f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val sphereObbDemo() {
    geomlib_obb box;
    box.center=geomlib_vec3_create(0,0,0);
    box.axes[0]=geomlib_vec3_create(1,0,0);
    box.axes[1]=geomlib_vec3_create(0,1,0);
    box.axes[2]=geomlib_vec3_create(0,0,1);
    box.half_extents=geomlib_vec3_create(1,1,1);
    geomlib_sphere_obb_result r;
    int hit=geomlib_sphere_obb_overlap(
        geomlib_vec3_create(2,0,0),1.5f,&box,&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.5f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val tetraTetraDemo() {
    int hit=geomlib_tetrahedra_intersect(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),geomlib_vec3_create(0,0,1),
        geomlib_vec3_create(0.5f,0,0),geomlib_vec3_create(1.5f,0,0),
        geomlib_vec3_create(0.5f,1,0),geomlib_vec3_create(0.5f,0,1));
    val result=val::object();
    result.set("success",hit==1);
    return result;
}

val sphereCapsuleDemo() {
    geomlib_sphere_capsule_result r;
    int hit=geomlib_sphere_capsule_overlap(
        geomlib_vec3_create(1,0,0),0.5f,
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(2,0,0),0.5f,&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+1.0f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val triTetraDemo() {
    int hit=geomlib_triangle_tetrahedron_intersect(
        geomlib_vec3_create(-1,0.3f,0.3f),geomlib_vec3_create(2,0.3f,0.3f),
        geomlib_vec3_create(-1,0.8f,0.3f),
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),geomlib_vec3_create(0,0,1));
    val result=val::object();
    result.set("success",hit==1);
    return result;
}

val polygonDistance2dDemo() {
    std::vector<geomlib_vec2> a={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)};
    std::vector<geomlib_vec2> b={
        geomlib_vec2_create(2,0),geomlib_vec2_create(3,0),
        geomlib_vec2_create(3,1),geomlib_vec2_create(2,1)};
    geomlib_polygon_distance2d_result r;
    int ok=geomlib_polygon_distance2d(a.data(),a.size(),b.data(),b.size(),&r);
    val result=val::object();
    result.set("success",ok && std::fabs(r.distance-1.0f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val capsuleAabbDemo() {
    geomlib_capsule_aabb_result r;
    int hit=geomlib_capsule_aabb_overlap(
        geomlib_vec3_create(2,0.5f,0.5f),geomlib_vec3_create(3,0.5f,0.5f),
        1.5f,geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,1,1),&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.5f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val polylineClosestDemo() {
    std::vector<geomlib_vec3> p={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(2,0,0),
        geomlib_vec3_create(2,2,0)};
    geomlib_polyline_closest r;
    int ok=geomlib_polyline_closest_point(
        p.data(),p.size(),geomlib_vec3_create(1,1,0),&r);
    val result=val::object();
    result.set("success",ok && std::fabs(r.distance-1.0f)<1e-3 &&
                         std::fabs(r.point.x-1.0f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val rayTetraDemo() {
    geomlib_ray_tetra_hit h;
    int hit=geomlib_ray_tetrahedron_intersect(
        geomlib_vec3_create(-1,-1,-1),geomlib_vec3_create(1,1,1),
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),geomlib_vec3_create(0,0,1),&h);
    val result=val::object();
    result.set("success",hit && std::fabs(h.t-1.0f)<1e-3 &&
                         std::fabs(h.point.x)<1e-3);
    result.set("t",static_cast<double>(h.t));
    return result;
}

val circleLineDemo() {
    geomlib_circle_line_intersection r=geomlib_circle_line_intersect(
        geomlib_vec2_create(-2,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(0,0),1.0f);
    val result=val::object();
    result.set("success",r.kind==GEOMLIB_CIRCLE_LINE_TWO_POINTS &&
                         std::fabs(r.p0.x+1.0f)<1e-3 &&
                         std::fabs(r.p1.x-1.0f)<1e-3);
    result.set("kind",static_cast<int>(r.kind));
    return result;
}

val planeSphereDemo() {
    geomlib_plane_sphere_intersection r=geomlib_plane_sphere_intersect(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),
        geomlib_vec3_create(0,0,1),2.0f);
    val result=val::object();
    result.set("success",r.kind==GEOMLIB_PLANE_SPHERE_CIRCLE &&
                         std::fabs(r.radius-1.7320508f)<1e-3);
    result.set("radius",static_cast<double>(r.radius));
    return result;
}

val triTriDistanceDemo() {
    geomlib_tri_tri_distance_result r;
    int ok=geomlib_triangle_triangle_distance(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),
        geomlib_vec3_create(3,0,0),geomlib_vec3_create(4,0,0),
        geomlib_vec3_create(3,1,0),&r);
    val result=val::object();
    result.set("success",ok && std::fabs(r.distance-2.0f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val rayObbDemo() {
    geomlib_obb box;
    box.center=geomlib_vec3_create(0,0,0);
    box.axes[0]=geomlib_vec3_create(1,0,0);
    box.axes[1]=geomlib_vec3_create(0,1,0);
    box.axes[2]=geomlib_vec3_create(0,0,1);
    box.half_extents=geomlib_vec3_create(1,1,1);
    geomlib_ray_obb_hit h;
    int hit=geomlib_ray_obb_intersect(
        geomlib_vec3_create(-2,0,0),geomlib_vec3_create(1,0,0),&box,&h);
    val result=val::object();
    result.set("success",hit && std::fabs(h.t-1.0)<1e-3 &&
                         std::fabs(h.normal.x+1.0f)<1e-3);
    result.set("t",h.t);
    return result;
}

val aabbDistanceDemo() {
    geomlib_aabb_distance_result r;
    int ok=geomlib_aabb_distance(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,1,1),
        geomlib_vec3_create(2,0,0),geomlib_vec3_create(3,1,1),&r);
    val result=val::object();
    result.set("success",ok && std::fabs(r.distance-1.0f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val polygonWidthDemo() {
    std::vector<geomlib_vec2> sq={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)};
    double w=geomlib_convex_polygon_min_width(sq.data(),sq.size());
    val result=val::object();
    result.set("success",std::fabs(w-1.0)<1e-4);
    result.set("width",w);
    return result;
}

val polygonDiameterDemo() {
    std::vector<geomlib_vec2> sq={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)};
    double d=geomlib_convex_polygon_diameter(sq.data(),sq.size());
    val result=val::object();
    result.set("success",std::fabs(d-1.41421356)<1e-4);
    result.set("diameter",d);
    return result;
}

val triAabbDemo() {
    int hit=geomlib_triangle_aabb_overlap(
        geomlib_vec3_create(0.5f,0.5f,-1),
        geomlib_vec3_create(0.5f,1.5f,0.5f),
        geomlib_vec3_create(1.5f,0.5f,0.5f),
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,1,1));
    val result=val::object();
    result.set("success",hit==1);
    return result;
}

val sphereTriangleDemo() {
    geomlib_sphere_triangle_result r;
    int hit=geomlib_sphere_triangle_overlap(
        geomlib_vec3_create(0.5f,0.5f,0.1f),0.2f,
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(2,0,0),
        geomlib_vec3_create(0,2,0),&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.1f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val capsuleDemo() {
    geomlib_capsule_result r;
    int hit=geomlib_capsule_capsule_overlap(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,2),0.5f,
        geomlib_vec3_create(0,0,3),geomlib_vec3_create(0,0,5),0.5f,&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val rayCylinderDemo() {
    geomlib_ray_cylinder_hit h;
    int hit=geomlib_ray_cylinder_intersect(
        geomlib_vec3_create(-2,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),1.0f,&h);
    val result=val::object();
    result.set("success",hit && std::fabs(h.t-1.0)<1e-3 &&
                         std::fabs(h.normal.x+1.0f)<1e-3);
    result.set("t",h.t);
    return result;
}

val segmentTriangleDemo() {
    geomlib_segment_triangle_hit h;
    int hit=geomlib_segment_triangle_intersect(
        geomlib_vec3_create(0.2f,0.2f,-1),geomlib_vec3_create(0.2f,0.2f,1),
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),&h);
    val result=val::object();
    result.set("success",hit && std::fabs(h.t-0.5f)<1e-3 &&
                         std::fabs(h.bary.x-0.6f)<1e-3);
    result.set("t",static_cast<double>(h.t));
    return result;
}

val circleCircleDemo() {
    geomlib_circle2d_intersection r=geomlib_circle2d_intersect(
        geomlib_vec2_create(0,0),1.0f,geomlib_vec2_create(1,0),1.0f);
    val result=val::object();
    result.set("success",r.kind==GEOMLIB_CIRCLE2D_TWO_POINTS &&
                         std::fabs(r.p0.x-0.5f)<1e-3 &&
                         std::fabs(std::fabs(r.p0.y)-0.866025f)<1e-3);
    result.set("kind",static_cast<int>(r.kind));
    return result;
}

val sphereAabbDemo() {
    geomlib_sphere_aabb_result r;
    int hit=geomlib_sphere_aabb_overlap(
        geomlib_vec3_create(2,0.5f,0.5f),1.5f,
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,1,1),&r);
    val result=val::object();
    result.set("success",hit && std::fabs(r.distance+0.5f)<1e-3);
    result.set("distance",static_cast<double>(r.distance));
    return result;
}

val planePlaneDemo() {
    geomlib_plane_line l;
    int hit=geomlib_plane_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),
        geomlib_vec3_create(1,0,0),geomlib_vec3_create(1,0,0),&l);
    val result=val::object();
    result.set("success",hit && std::fabs(l.point.x-1.0f)<1e-3 &&
                         std::fabs(l.direction.y-1.0f)<1e-3);
    return result;
}

val raySphereDemo() {
    geomlib_ray_sphere_hit h;
    int hit=geomlib_ray_sphere_intersect(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(5,0,0),1.0f,&h);
    val result=val::object();
    result.set("success",hit && std::fabs(h.t0-4.0)<1e-4 &&
                         std::fabs(h.t1-6.0)<1e-4);
    result.set("t0",h.t0);
    return result;
}

val bezierDemo() {
    std::vector<geomlib_vec3> cp={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(2,0,0),geomlib_vec3_create(3,0,0)
    };
    geomlib_vec3 p=geomlib_bezier_evaluate(cp.data(),cp.size(),0.5);
    val result=val::object();
    result.set("success",std::fabs(p.x-1.5f)<1e-4);
    result.set("value",static_cast<double>(p.x));
    return result;
}

val vwDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i=0;i<8;++i)
        pts.push_back(geomlib_vec3_create((float)i,(i==4)?5.0f:0.0f,0));
    std::vector<uint32_t> idx(pts.size());
    size_t kept=0;
    int rc=geomlib_simplify_polyline_vw(pts.data(),pts.size(),0.1,
                                        idx.data(),&kept);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && kept==5 && idx[2]==4);
    result.set("kept",static_cast<int>(kept));
    return result;
}

val solidAngleDemo() {
    geomlib_vec3 o=geomlib_vec3_create(0,0,0);
    geomlib_vec3 ax=geomlib_vec3_create(1,0,0);
    geomlib_vec3 ay=geomlib_vec3_create(0,1,0);
    geomlib_vec3 az=geomlib_vec3_create(0,0,1);
    double omega=geomlib_solid_angle_triangle(o,ax,ay,az);
    geomlib_trimesh box;
    geomlib_create_box(&box,2.0f,2.0f,2.0f);
    double wn=geomlib_winding_number_mesh(&box,o);
    val result=val::object();
    result.set("success",std::fabs(omega-1.57079632679)<1e-4 &&
                         std::fabs(wn-1.0)<1e-3);
    result.set("omega",omega);
    geomlib_trimesh_free(&box);
    return result;
}

val meshOverlapDemo() {
    geomlib_trimesh a, b;
    geomlib_create_box(&a, 1.0f, 1.0f, 1.0f);
    geomlib_create_box(&b, 1.0f, 1.0f, 1.0f);
    for (size_t i = 0; i < b.vertex_count; ++i) b.vertices[i].x += 5.0f;
    geomlib_overlap_pairs p = geomlib_overlap_pairs_create();
    int rc = geomlib_aabb_overlap_pairs(&a, &b, &p);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && p.count == 0);
    result.set("pairs", static_cast<int>(p.count));
    geomlib_overlap_pairs_free(&p);
    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    return result;
}

val visibilityDemo() {
    geomlib_vec2 s=geomlib_vec2_create(-3,0),g=geomlib_vec2_create(3,0);
    geomlib_vec2 box[4]={
        geomlib_vec2_create(-1,-1),geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1),geomlib_vec2_create(-1,1)
    };
    const geomlib_vec2 *obs[1]={box};
    size_t sizes[1]={4};
    geomlib_visibility_path p=geomlib_visibility_path_create();
    int rc=geomlib_visibility_shortest_path(&s,&g,obs,sizes,1,&p);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && p.length>6.0f);
    result.set("length",static_cast<double>(p.length));
    geomlib_visibility_path_free(&p);
    return result;
}

val raster2dDemo() {
    std::vector<geomlib_vec2> sq={
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2),geomlib_vec2_create(0,2)
    };
    geomlib_raster2d r=geomlib_raster2d_create();
    int rc=geomlib_rasterize_polygon2d(sq.data(),sq.size(),30,30,0.1f,0,&r);
    int inside=0;
    if (rc==GEOMLIB_MESH_SUCCESS)
        for (size_t i=0;i<(size_t)r.nx*r.ny;++i) if (r.mask[i]) inside++;
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && inside>100);
    result.set("inside",inside);
    geomlib_raster2d_free(&r);
    return result;
}

val cage2dHarmonicDemo() {
    std::vector<geomlib_vec2> cage={
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2),geomlib_vec2_create(0,2)
    };
    std::vector<geomlib_vec2> moved={
        geomlib_vec2_create(1,-1),geomlib_vec2_create(3,-1),
        geomlib_vec2_create(3,1),geomlib_vec2_create(1,1)
    };
    geomlib_vec2 q=geomlib_vec2_create(0.7f,0.4f),out;
    int rc=geomlib_cage2d_harmonic_deform(cage.data(),moved.data(),4,&q,1,&out);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(out.x-(q.x+1.0f))<1e-3 &&
                         std::fabs(out.y-(q.y-1.0f))<1e-3);
    return result;
}

val curveSmoothDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i=0;i<11;++i)
        pts.push_back(geomlib_vec3_create((float)i,(i%2)?0.5f:0.0f,0));
    int rc=geomlib_smooth_polyline(pts.data(),pts.size(),30,0.5f);
    double rough=0.0;
    for (size_t i=1;i+1<pts.size();++i) {
        geomlib_vec3 d=geomlib_vec3_sub(pts[i],
            geomlib_vec3_mul(geomlib_vec3_add(pts[i-1],pts[i+1]),0.5f));
        rough+=geomlib_vec3_length_squared(d);
    }
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && rough<0.01);
    result.set("roughness",rough);
    return result;
}

val poisson2dDemo() {
    std::vector<geomlib_vec2> box={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)
    };
    geomlib_poisson2d_params p;
    p.radius=0.15f; p.candidates=30; p.seed=5; p.max_samples=200;
    std::vector<geomlib_vec2> pts(p.max_samples);
    size_t n=geomlib_poisson_disk_polygon(box.data(),box.size(),&p,
                                          pts.data(),pts.size());
    val result=val::object();
    result.set("success",n>10);
    result.set("count",static_cast<int>(n));
    return result;
}

val frechetDemo() {
    std::vector<geomlib_vec3> a={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0)};
    std::vector<geomlib_vec3> b={
        geomlib_vec3_create(0,2,0),geomlib_vec3_create(1,2,0)};
    double d=geomlib_frechet_distance(a.data(),a.size(),b.data(),b.size());
    val result=val::object();
    result.set("success",std::fabs(d-2.0)<1e-4);
    result.set("distance",d);
    return result;
}

val sphereIntersectionDemo() {
    geomlib_sphere_intersection r=geomlib_sphere_sphere_intersection(
        geomlib_vec3_create(0,0,0),2.0f,geomlib_vec3_create(3,0,0),2.0f);
    val result=val::object();
    result.set("success",r.kind==GEOMLIB_SPHERE_CIRCLE &&
                         std::fabs(r.center.x-1.5)<1e-3);
    result.set("radius",static_cast<double>(r.radius));
    return result;
}

val pcaDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i=0;i<100;++i) {
        float u=(float)((i*37)%101)/10.0f-5.0f;
        float v=(float)((i*53)%97)/100.0f-0.5f;
        pts.push_back(geomlib_vec3_create(u*3.0f,v,0));
    }
    geomlib_pca_result p;
    int rc=geomlib_compute_pca(pts.data(),pts.size(),&p);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(p.eigenvectors[0].x)>0.99);
    result.set("eig0",static_cast<double>(p.eigenvalues[0]));
    return result;
}

val mstDemo() {
    double pts[8]={0,0, 1,0, 1,1, 0,1};
    geomlib_mst_result r=geomlib_mst_result_create();
    int rc=geomlib_compute_mst(pts,4,2,&r);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         r.edge_count==3 && std::fabs(r.total_weight-3.0)<1e-4);
    result.set("weight",static_cast<double>(r.total_weight));
    geomlib_mst_result_free(&r);
    return result;
}

val segment2dDemo() {
    geomlib_seg2d_result r=geomlib_segment2d_intersect(
        geomlib_vec2_create(-1,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(0,-1),geomlib_vec2_create(0,1));
    val result=val::object();
    result.set("success",r.kind==GEOMLIB_SEG2D_PROPER && std::fabs(r.point.x)<1e-4);
    result.set("x",static_cast<double>(r.point.x));
    return result;
}

val meshClipDemo() {
    geomlib_trimesh box;
    int rc=geomlib_create_box(&box,1.0f,1.0f,1.0f);
    geomlib_trimesh out;
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_clip_mesh_plane(&box,geomlib_vec3_create(0,0,0),
            geomlib_vec3_create(1,0,0),0,&out);
    geomlib_mass_properties mp;
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_compute_mass_properties(&out,&mp);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(mp.volume-0.5)<0.02);
    result.set("volume",static_cast<double>(mp.volume));
    if (out.triangle_count) geomlib_trimesh_free(&out);
    geomlib_trimesh_free(&box);
    return result;
}

val quaternionDemo() {
    geomlib_quat a=geomlib_quat_from_axis_angle(
        geomlib_vec3_create(0,0,1),1.0);
    geomlib_quat mid=geomlib_quat_slerp(geomlib_quat_identity(),a,0.5);
    geomlib_vec3 v=geomlib_quat_rotate(mid,geomlib_vec3_create(1,0,0));
    val result=val::object();
    result.set("success",std::fabs(geomlib_vec3_length(v)-1.0)<1e-4 &&
                         std::fabs(v.y)>0.4);
    result.set("y",static_cast<double>(v.y));
    return result;
}

val laplacianDemo() {
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid,9,8);
    for (int y=0;y<3;++y)
        for (int x=0;x<3;++x)
            geomlib_trimesh_add_vertex(&grid,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<2;++y)
        for (int x=0;x<2;++x) {
            uint32_t i=(uint32_t)(y*3+x);
            geomlib_trimesh_add_triangle(&grid,i,i+1,i+4);
            geomlib_trimesh_add_triangle(&grid,i,i+4,i+3);
        }
    geomlib_laplacian l=geomlib_laplacian_create();
    int rc=geomlib_build_laplacian(&grid,&l);
    std::vector<double> ones(l.vertex_count,1.0),y(l.vertex_count,0.0);
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_laplacian_apply(&l,ones.data(),0.0,y.data());
    double mx=0.0;
    for (double v:y) if (std::fabs(v)>mx) mx=std::fabs(v);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && mx<1e-9);
    result.set("max",mx);
    geomlib_laplacian_free(&l);
    geomlib_trimesh_free(&grid);
    return result;
}

val pointInPolygonDemo() {
    std::vector<geomlib_vec2> sq={
        geomlib_vec2_create(-1,-1),geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1),geomlib_vec2_create(-1,1)
    };
    int inside=geomlib_point_in_polygon(sq.data(),sq.size(),
                                        geomlib_vec2_create(0.2f,0.3f));
    int outside=geomlib_point_in_polygon(sq.data(),sq.size(),
                                         geomlib_vec2_create(2,2));
    val result=val::object();
    result.set("success",inside==1 && outside==0);
    return result;
}

val poisson3dDemo() {
    geomlib_poisson3d_params p;
    p.box_min=geomlib_vec3_create(0,0,0);
    p.box_max=geomlib_vec3_create(1,1,1);
    p.radius=0.2f;
    p.candidates=20;
    p.seed=42;
    p.max_samples=500;
    std::vector<geomlib_vec3> pts(p.max_samples);
    size_t n=geomlib_poisson_disk_volume(&p,pts.data(),pts.size());
    val result=val::object();
    result.set("success",n>10);
    result.set("count",static_cast<int>(n));
    return result;
}

val cage3dDemo() {
    geomlib_vec3 cage[4]={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),geomlib_vec3_create(0,0,1)
    };
    geomlib_vec3 moved[4];
    for (int i=0;i<4;++i)
        moved[i]=geomlib_vec3_add(geomlib_vec3_mul(cage[i],2.0f),
                                  geomlib_vec3_create(1,2,3));
    geomlib_vec3 q=geomlib_vec3_create(0.2f,0.3f,0.1f),out;
    int rc=geomlib_deform_tetra_cage(cage,moved,&q,1,&out);
    geomlib_vec3 expect=geomlib_vec3_add(geomlib_vec3_mul(q,2.0f),
                                        geomlib_vec3_create(1,2,3));
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         geomlib_vec3_length(geomlib_vec3_sub(out,expect))<1e-4f);
    return result;
}

val geodesicDiameterDemo() {
    geomlib_trimesh m;
    int rc=geomlib_create_sphere(&m,1.0f,10,7);
    geomlib_geodesic_diameter d;
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_geodesic_diameter_sweep(&m,&d);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && d.distance>2.5f);
    result.set("distance",static_cast<double>(d.distance));
    geomlib_trimesh_free(&m);
    return result;
}

val polyboolHolesDemo() {
    std::vector<geomlib_vec2> outer={
        geomlib_vec2_create(0,0),geomlib_vec2_create(4,0),
        geomlib_vec2_create(4,4),geomlib_vec2_create(0,4)
    };
    std::vector<geomlib_vec2> hole={
        geomlib_vec2_create(1,1),geomlib_vec2_create(3,1),
        geomlib_vec2_create(3,3),geomlib_vec2_create(1,3)
    };
    const geomlib_vec2 *holes[1]={hole.data()};
    size_t sizes[1]={hole.size()};
    geomlib_polygon_with_holes A,B;
    int rc=geomlib_polygon_with_holes_set(outer.data(),outer.size(),holes,sizes,1,&A);
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_polygon_with_holes_set(outer.data(),outer.size(),nullptr,nullptr,0,&B);
    geomlib_poly_bool_result r=geomlib_poly_bool_result_create();
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_polygon_boolean_holes(&A,&B,GEOMLIB_POLY_INTERSECTION,&r);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS);
    result.set("loops",static_cast<int>(r.loop_count));
    geomlib_poly_bool_result_free(&r);
    geomlib_polygon_with_holes_free(&A);
    geomlib_polygon_with_holes_free(&B);
    return result;
}

val triTriDemo() {
    geomlib_tri_tri_result r=geomlib_triangle_intersection(
        geomlib_vec3_create(-1,-1,0),geomlib_vec3_create(1,-1,0),
        geomlib_vec3_create(0,1,0),
        geomlib_vec3_create(-1,1,0),geomlib_vec3_create(1,1,0),
        geomlib_vec3_create(0,-1,0));
    val result=val::object();
    result.set("success",r.intersects==1);
    return result;
}

val rayPlaneDemo() {
    float t;
    geomlib_vec3 hit;
    int ok=geomlib_ray_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),
        geomlib_vec3_create(0,0,3),geomlib_vec3_create(0,0,1),&t,&hit);
    val result=val::object();
    result.set("success",ok && std::fabs(t-3.0)<1e-4);
    result.set("t",static_cast<double>(t));
    return result;
}

val segmentDistanceDemo() {
    geomlib_segment_distance d;
    int rc=geomlib_segment_closest_points(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(4,0,0),
        geomlib_vec3_create(2,3,0),geomlib_vec3_create(2,5,0),&d);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(d.distance-3.0)<1e-4);
    result.set("distance",static_cast<double>(d.distance));
    return result;
}

val cylinderFitDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i=0;i<120;++i) {
        float t=(float)((i%20)-10)*0.3f;
        float a=(float)(i/20)*6.283185f/6.0f;
        pts.push_back(geomlib_vec3_create(1.0f+2.0f*std::cos(a),
                                          2.0f+2.0f*std::sin(a),t));
    }
    geomlib_cylinder_fit fit;
    int rc=geomlib_fit_cylinder(pts.data(),pts.size(),
                                geomlib_vec3_create(0,0,1),&fit);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(fit.radius-2.0)<0.05);
    result.set("radius",static_cast<double>(fit.radius));
    return result;
}

val ellipseFitDemo() {
    std::vector<geomlib_vec2> pts;
    for (int i=0;i<80;++i) {
        float t=(float)(i*6.283185/80.0);
        float u=2.0f*std::cos(t),v=1.0f*std::sin(t);
        pts.push_back(geomlib_vec2_create(u,v));
    }
    geomlib_ellipse_fit fit;
    int rc=geomlib_fit_ellipse(pts.data(),pts.size(),&fit);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(fit.axis_a-2.0)<0.02 &&
                         std::fabs(fit.axis_b-1.0)<0.02);
    result.set("a",static_cast<double>(fit.axis_a));
    return result;
}

val chaikinDemo() {
    std::vector<geomlib_vec3> line={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(2,0,0)
    };
    std::vector<geomlib_vec3> out(5);
    size_t n=geomlib_chaikin_subdivide(line.data(),line.size(),1,out.data(),out.size());
    val result=val::object();
    result.set("success",n==5 && out[0].x==0.0f && out[4].x==2.0f);
    result.set("count",static_cast<int>(n));
    return result;
}

val planeFitDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i=0;i<50;++i)
        pts.push_back(geomlib_vec3_create((float)i,(float)(i%7),2.0f));
    geomlib_plane_fit fit;
    int rc=geomlib_fit_plane(pts.data(),pts.size(),&fit);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(fit.normal.z)>0.99 && fit.rms<1e-3);
    result.set("rms",static_cast<double>(fit.rms));
    return result;
}

val sphereFitDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i=0;i<120;++i) {
        float u=(float)(i*2.399963);
        float z=1.0f-2.0f*(float)(i+0.5)/120.0f;
        float r=std::sqrt(std::fmax(0.0f,1.0f-z*z));
        pts.push_back(geomlib_vec3_create(r*std::cos(u),r*std::sin(u),z));
    }
    geomlib_sphere_fit fit;
    int rc=geomlib_fit_sphere(pts.data(),pts.size(),&fit);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(fit.radius-1.0)<1e-2);
    result.set("radius",static_cast<double>(fit.radius));
    return result;
}

val circleFitDemo() {
    std::vector<geomlib_vec2> pts;
    for (int i=0;i<60;++i) {
        float a=(float)(i*6.283185/60.0);
        pts.push_back(geomlib_vec2_create(1.0f+2.0f*std::cos(a),
                                          2.0f+2.0f*std::sin(a)));
    }
    geomlib_circle_fit fit;
    int rc=geomlib_fit_circle(pts.data(),pts.size(),&fit);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(fit.radius-2.0)<1e-3);
    result.set("radius",static_cast<double>(fit.radius));
    return result;
}

val rigidFitDemo() {
    std::vector<geomlib_vec3> src,dst;
    for (int i=0;i<20;++i) {
        src.push_back(geomlib_vec3_create((float)i,(float)(i%5),0));
        dst.push_back(geomlib_vec3_create((float)i+1.5f,(float)(i%5)-1.0f,0));
    }
    geomlib_mat4 fit;
    float rms=-1.0f;
    int rc=geomlib_fit_rigid_points(src.data(),dst.data(),src.size(),false,&fit,&rms);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && rms<1e-3);
    result.set("rms",static_cast<double>(rms));
    return result;
}

val geodesicBallDemo() {
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&m,i,i+6,i+5);
        }
    geomlib_geodesic_ball b=geomlib_geodesic_ball_create();
    int rc=geomlib_extract_geodesic_ball(&m,0,2.0f,&b);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && b.count>1 && b.count<25);
    result.set("count",static_cast<int>(b.count));
    geomlib_geodesic_ball_free(&b);
    geomlib_trimesh_free(&m);
    return result;
}

val geodesicPathDemo() {
    geomlib_trimesh m;
    int rc=geomlib_create_sphere(&m,1.0f,10,7);
    geomlib_geodesic_path p=geomlib_geodesic_path_create();
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_geodesic_path_between(&m,0,(uint32_t)(m.vertex_count/2),&p);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && p.count>1);
    result.set("vertices",static_cast<int>(p.count));
    result.set("length",static_cast<double>(p.length));
    geomlib_geodesic_path_free(&p);
    geomlib_trimesh_free(&m);
    return result;
}

val minCircleDemo() {
    std::vector<geomlib_vec2> pts;
    for (int i=0;i<60;++i) {
        float a=(float)(i*6.283185/60.0);
        pts.push_back(geomlib_vec2_create(std::cos(a),std::sin(a)));
    }
    geomlib_min_circle2d c;
    int rc=geomlib_minimum_enclosing_circle(pts.data(),pts.size(),&c);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(c.radius-1.0)<1e-3);
    result.set("radius",static_cast<double>(c.radius));
    return result;
}

val polygonClipDemo() {
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1)
    };
    std::vector<geomlib_vec2> half = {
        geomlib_vec2_create(-10,-10), geomlib_vec2_create(0,-10),
        geomlib_vec2_create(0,10), geomlib_vec2_create(-10,10)
    };
    geomlib_polygon2d out={};
    int rc=geomlib_clip_polygon_convex(sq.data(),sq.size(),half.data(),half.size(),&out);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && out.vertex_count==4);
    result.set("vertices",static_cast<int>(out.vertex_count));
    geomlib_polygon2d_free(&out);
    return result;
}

val polygonMomentsDemo() {
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1)
    };
    geomlib_polygon_moments m;
    int rc=geomlib_compute_polygon_moments(sq.data(),sq.size(),&m);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(m.area-1.0)<1e-6 &&
                         std::fabs(m.ixx-1.0/12.0)<1e-4);
    result.set("area",m.area);
    return result;
}

val resampleDemo() {
    std::vector<geomlib_vec3> line={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(2,0,0),
        geomlib_vec3_create(5,0,0),geomlib_vec3_create(10,0,0)
    };
    std::vector<geomlib_vec3> out(11);
    int rc=geomlib_resample_polyline(line.data(),line.size(),out.size(),out.data());
    double len=0.0;
    for (size_t i=1;i<out.size();++i)
        len+=geomlib_vec3_length(geomlib_vec3_sub(out[i],out[i-1]));
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(len-10.0)<1e-3);
    result.set("length",len);
    return result;
}

val dbscanDemo() {
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud,60);
    for (int i=0;i<30;++i)
        geomlib_pointcloud_add_point(&cloud,
            geomlib_vec3_create(-2.0f+0.05f*sinf((float)i),0,0),nullptr,1.0f);
    for (int i=0;i<30;++i)
        geomlib_pointcloud_add_point(&cloud,
            geomlib_vec3_create(2.0f+0.05f*sinf((float)i),0,0),nullptr,1.0f);
    geomlib_dbscan_params p=geomlib_dbscan_params_default();
    p.eps=0.5f; p.min_points=5;
    std::vector<uint32_t> labels(cloud.point_count);
    size_t clusters=0;
    int rc=geomlib_pointcloud_dbscan(&cloud,&p,labels.data(),&clusters);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && clusters==2);
    result.set("clusters",static_cast<int>(clusters));
    geomlib_pointcloud_free(&cloud);
    return result;
}

val kmeansDemo() {
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud,100);
    for (int i=0;i<100;++i) {
        float x=(i<50)?-2.0f:2.0f;
        geomlib_pointcloud_add_point(&cloud,geomlib_vec3_create(x,0,0),nullptr,1.0f);
    }
    std::vector<uint32_t> labels(cloud.point_count);
    int k=geomlib_pointcloud_kmeans(&cloud,2,20,labels.data(),nullptr);
    int c0=0,c1=0;
    for (uint32_t l:labels) { if(l==0)c0++; else if(l==1)c1++; }
    val result=val::object();
    result.set("success",k==2 && c0==50 && c1==50);
    result.set("clusters",static_cast<int>(c0+c1));
    geomlib_pointcloud_free(&cloud);
    return result;
}

val symmetryDemo() {
    geomlib_trimesh box;
    int rc=geomlib_create_box(&box,2.0f,1.0f,0.5f);
    geomlib_symmetry_result s;
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_detect_symmetry(&box,&s);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && s.rms_error<0.02f);
    result.set("rms",static_cast<double>(s.rms_error));
    geomlib_trimesh_free(&box);
    return result;
}

val meshDistanceDemo() {
    geomlib_trimesh a,b;
    int rc=geomlib_create_sphere(&a,1.0f,10,7);
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_create_sphere(&b,1.0f,10,7);
    if (rc==GEOMLIB_MESH_SUCCESS)
        for (size_t i=0;i<b.vertex_count;++i) b.vertices[i].x+=3.0f;
    geomlib_mesh_distance d;
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_mesh_mesh_distance(&a,&b,&d);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(d.distance-1.0)<0.06);
    result.set("distance",static_cast<double>(d.distance));
    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    return result;
}

val spectralSegmentDemo() {
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid,49,72);
    for (int y=0;y<7;++y)
        for (int x=0;x<7;++x)
            geomlib_trimesh_add_vertex(&grid,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<6;++y)
        for (int x=0;x<6;++x) {
            uint32_t i=(uint32_t)(y*7+x);
            geomlib_trimesh_add_triangle(&grid,i,i+1,i+8);
            geomlib_trimesh_add_triangle(&grid,i,i+8,i+7);
        }
    std::vector<uint32_t> labels(grid.vertex_count,0);
    int rc=geomlib_spectral_segment(&grid,2,4,labels.data());
    int c0=0,c1=0;
    for (uint32_t l:labels) { if(l==0)c0++; else if(l==1)c1++; }
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && c0>0 && c1>0);
    result.set("clusters",c0+c1);
    geomlib_trimesh_free(&grid);
    return result;
}

val ambientOcclusionDemo() {
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid,25,32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&grid,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&grid,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&grid,i,i+6,i+5);
        }
    std::vector<float> ao(grid.vertex_count);
    int rc=geomlib_mesh_ambient_occlusion(&grid,16,10.0f,ao.data());
    float mn=1.0f;
    for (float v:ao) if (v<mn) mn=v;
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && mn>0.95f);
    result.set("min",static_cast<double>(mn));
    geomlib_trimesh_free(&grid);
    return result;
}

val polylineIntersectionDemo() {
    geomlib_trimesh sphere;
    int rc=geomlib_create_sphere(&sphere,1.0f,10,7);
    geomlib_vec3 line[2]={geomlib_vec3_create(-2,0,0),geomlib_vec3_create(2,0,0)};
    geomlib_polyline_hits h=geomlib_polyline_hits_create();
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_polyline_mesh_intersection(line,2,&sphere,&h);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && h.count==2);
    result.set("hits",static_cast<int>(h.count));
    geomlib_polyline_hits_free(&h);
    geomlib_trimesh_free(&sphere);
    return result;
}

val minAreaRectDemo() {
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1)
    };
    geomlib_min_area_rect r;
    int rc=geomlib_minimum_area_rectangle(sq.data(),sq.size(),&r);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(r.area-4.0)<1e-3);
    result.set("area",static_cast<double>(r.area));
    return result;
}

val inscribedCircleDemo() {
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1)
    };
    geomlib_inscribed_circle c;
    int rc=geomlib_largest_inscribed_circle(sq.data(),sq.size(),0.08f,&c);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && std::fabs(c.radius-1.0f)<0.02);
    result.set("radius",static_cast<double>(c.radius));
    return result;
}

val massPropertiesDemo() {
    geomlib_trimesh box;
    int rc=geomlib_create_box(&box,1.0f,1.0f,1.0f);
    geomlib_mass_properties mp;
    if (rc==GEOMLIB_MESH_SUCCESS) rc=geomlib_compute_mass_properties(&box,&mp);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(mp.volume-1.0)<1e-3 &&
                         std::fabs(mp.inertia.m[0]-1.0/6.0)<0.01);
    result.set("volume",static_cast<double>(mp.volume));
    result.set("inertia",static_cast<double>(mp.inertia.m[0]));
    geomlib_trimesh_free(&box);
    return result;
}

val miniballDemo() {
    std::vector<geomlib_vec3> pts;
    for (int z=-1;z<=1;z+=2)
        for (int y=-1;y<=1;y+=2)
            for (int x=-1;x<=1;x+=2)
                pts.push_back(geomlib_vec3_create((float)x,(float)y,(float)z));
    geomlib_bounding_sphere s;
    int rc=geomlib_minimum_enclosing_sphere(pts.data(),pts.size(),&s);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(s.radius-1.73205)<1e-3);
    result.set("radius",static_cast<double>(s.radius));
    return result;
}

val boundingSphereDemo() {
    std::vector<geomlib_vec3> pts;
    for (int z=-1;z<=1;z+=2)
        for (int y=-1;y<=1;y+=2)
            for (int x=-1;x<=1;x+=2)
                pts.push_back(geomlib_vec3_create((float)x,(float)y,(float)z));
    geomlib_bounding_sphere s;
    int rc=geomlib_bounding_sphere_points(pts.data(),pts.size(),&s);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS &&
                         std::fabs(s.radius-1.73205)<0.05);
    result.set("radius",static_cast<double>(s.radius));
    return result;
}

val convexDecomposeDemo() {
    std::vector<geomlib_vec2> L = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,1), geomlib_vec2_create(1,1),
        geomlib_vec2_create(1,2), geomlib_vec2_create(0,2)
    };
    geomlib_poly_decomposition d=geomlib_poly_decomposition_create();
    int rc=geomlib_convex_decompose2d(L.data(),L.size(),&d);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && d.part_count>=2);
    result.set("parts",static_cast<int>(d.part_count));
    geomlib_poly_decomposition_free(&d);
    return result;
}

val hksDemo() {
    geomlib_trimesh sphere;
    int rc=geomlib_create_sphere(&sphere,1.0f,10,7);
    double t=0.5;
    std::vector<double> hks(sphere.vertex_count,0.0);
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_compute_hks(&sphere,8,&t,1,hks.data());
    double mx=0.0;
    for (double v:hks) if (v>mx) mx=v;
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && mx>0.0 && mx<=1.000001);
    result.set("count",static_cast<int>(hks.size()));
    result.set("max",mx);
    geomlib_trimesh_free(&sphere);
    return result;
}

val pointCurvatureDemo() {
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud,200);
    for (int i=0;i<200;++i)
        geomlib_pointcloud_add_point(&cloud,
            geomlib_vec3_create(((i*37)%101)/100.0f,((i*73)%97)/100.0f,0),
            nullptr,1.0f);
    std::vector<float> curv(cloud.point_count);
    int rc=geomlib_pointcloud_curvature(&cloud,16,curv.data(),nullptr);
    float mx=0.0f;
    for (float c:curv) if (c>mx) mx=c;
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && mx<0.01f);
    result.set("max",static_cast<double>(mx));
    geomlib_pointcloud_free(&cloud);
    return result;
}

val regionGrowingDemo() {
    geomlib_trimesh box;
    int rc=geomlib_create_box(&box,1.0f,1.0f,1.0f);
    std::vector<uint32_t> labels(box.triangle_count);
    size_t regions=0;
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_mesh_region_growing(&box,30.0f,labels.data(),&regions);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && regions==6);
    result.set("regions",static_cast<int>(regions));
    geomlib_trimesh_free(&box);
    return result;
}

val meshContourDemo() {
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&m,i,i+6,i+5);
        }
    std::vector<float> values(m.vertex_count);
    for (size_t i=0;i<m.vertex_count;++i) values[i]=m.vertices[i].x;
    geomlib_mesh_contour c=geomlib_mesh_contour_create();
    int rc=geomlib_mesh_contour_lines(&m,values.data(),2.0f,&c);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && c.segment_count>0);
    result.set("segments",static_cast<int>(c.segment_count));
    geomlib_mesh_contour_free(&c);
    geomlib_trimesh_free(&m);
    return result;
}

val icpPlaneDemo() {
    geomlib_trimesh target, source;
    int rc=geomlib_create_sphere(&target,1.0f,10,7);
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_trimesh_init(&source,target.vertex_count,target.triangle_count);
    if (rc==GEOMLIB_MESH_SUCCESS) {
        for (size_t i=0;i<target.vertex_count;++i)
            geomlib_trimesh_add_vertex(&source,
                geomlib_vec3_create(target.vertices[i].x+0.15f,
                                    target.vertices[i].y+0.1f,
                                    target.vertices[i].z),nullptr,nullptr);
        for (size_t i=0;i<target.triangle_count;++i)
            geomlib_trimesh_add_triangle(&source,target.indices[i*3],
                target.indices[i*3+1],target.indices[i*3+2]);
    }
    geomlib_mat4 T=geomlib_mat4_identity();
    float rms=1.0f;
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_icp_point_to_plane(&source,&target,nullptr,&T,&rms);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && rms<0.01f);
    result.set("rms",static_cast<double>(rms));
    if (source.vertex_count) geomlib_trimesh_free(&source);
    if (target.vertex_count) geomlib_trimesh_free(&target);
    return result;
}

val surfaceVoronoiDemo() {
    geomlib_trimesh m;
    int rc=geomlib_create_sphere(&m,1.0f,10,7);
    uint32_t sites[2]={0,30};
    geomlib_surface_voronoi v=geomlib_surface_voronoi_create();
    if (rc==GEOMLIB_MESH_SUCCESS && sites[1]>=m.vertex_count)
        sites[1]=(uint32_t)(m.vertex_count/2);
    if (rc==GEOMLIB_MESH_SUCCESS)
        rc=geomlib_surface_voronoi_compute(&m,sites,2,&v);
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && v.sizes[0]>0 && v.sizes[1]>0);
    result.set("first",static_cast<int>(v.sizes[0]));
    result.set("second",static_cast<int>(v.sizes[1]));
    geomlib_surface_voronoi_free(&v);
    if (m.vertex_count) geomlib_trimesh_free(&m);
    return result;
}

val refine2dDemo() {
    std::vector<geomlib_vec2> rect = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(3,0),
        geomlib_vec2_create(3,1), geomlib_vec2_create(0,1)
    };
    geomlib_refine2d_params p = geomlib_refine2d_params_default();
    p.min_angle_deg = 25.0f;
    p.max_points = 40;
    geomlib_cdt2d out = geomlib_cdt2d_create();
    int rc = geomlib_refine_polygon2d(rect.data(), rect.size(), &p, &out);
    double area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < out.triangle_count; ++i) {
            geomlib_vec2 a=out.vertices[out.triangles[i*3]];
            geomlib_vec2 b=out.vertices[out.triangles[i*3+1]];
            geomlib_vec2 c=out.vertices[out.triangles[i*3+2]];
            area += 0.5*std::fabs((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
        }
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS &&
                          out.vertex_count > 4 && std::fabs(area-3.0) < 1e-3);
    result.set("points", static_cast<int>(out.vertex_count));
    geomlib_cdt2d_free(&out);
    return result;
}

val collapseDemo() {
    geomlib_trimesh m;
    int rc=geomlib_create_sphere(&m,1.0f,8,6);
    if (rc==GEOMLIB_MESH_SUCCESS) geomlib_refine_edges(&m,0.45f,150);
    size_t before=m.triangle_count;
    int collapses = rc==GEOMLIB_MESH_SUCCESS
                  ? geomlib_coarsen_edges(&m,0.25f,200,1) : 0;
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && collapses>0 &&
                         m.triangle_count<before);
    result.set("collapses",collapses);
    result.set("triangles",static_cast<int>(m.triangle_count));
    geomlib_trimesh_free(&m);
    return result;
}

val refineDemo() {
    geomlib_trimesh m;
    int rc=geomlib_create_sphere(&m,1.0f,8,6);
    int splits=0;
    float before=0.0f,after=0.0f;
    if (rc==GEOMLIB_MESH_SUCCESS) {
        for (size_t i=0;i<m.triangle_count;++i)
            for (int k=0;k<3;++k) {
                uint32_t a=m.indices[i*3+k],b=m.indices[i*3+(k+1)%3];
                before=std::fmax(before,geomlib_vec3_length(
                    geomlib_vec3_sub(m.vertices[a],m.vertices[b])));
            }
        splits=geomlib_refine_edges(&m,0.45f,150);
        for (size_t i=0;i<m.triangle_count;++i)
            for (int k=0;k<3;++k) {
                uint32_t a=m.indices[i*3+k],b=m.indices[i*3+(k+1)%3];
                after=std::fmax(after,geomlib_vec3_length(
                    geomlib_vec3_sub(m.vertices[a],m.vertices[b])));
            }
    }
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && splits>0 && after<before);
    result.set("splits",splits);
    result.set("after",static_cast<double>(after));
    geomlib_trimesh_free(&m);
    return result;
}

val qualityFlipDemo() {
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x) {
            float xx=(float)x+((x>0&&x<4&&y>0&&y<4)?0.45f*sinf((float)(x*3+y)):0.0f);
            float yy=(float)y+((x>0&&x<4&&y>0&&y<4)?0.30f*cosf((float)(x*5+y*2)):0.0f);
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create(xx,yy,0),nullptr,nullptr);
        }
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            if ((x+y)%2) { geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
                           geomlib_trimesh_add_triangle(&m,i,i+6,i+5); }
            else { geomlib_trimesh_add_triangle(&m,i,i+1,i+5);
                   geomlib_trimesh_add_triangle(&m,i+1,i+6,i+5); }
        }
    float before=geomlib_mesh_min_angle(&m);
    int flips=geomlib_improve_mesh_quality(&m,30);
    float after=geomlib_mesh_min_angle(&m);
    val result=val::object();
    result.set("success",flips>0 && after>before);
    result.set("flips",flips);
    result.set("before",static_cast<double>(before));
    result.set("after",static_cast<double>(after));
    geomlib_trimesh_free(&m);
    return result;
}

val splineDemo() {
    std::vector<geomlib_vec3> pts = {
        geomlib_vec3_create(0,0,0), geomlib_vec3_create(1,0.3f,0),
        geomlib_vec3_create(2,-0.2f,0), geomlib_vec3_create(3,0,0)
    };
    size_t count=geomlib_catmull_output_count(pts.size(),0,8);
    std::vector<geomlib_vec3> out(count);
    int rc=geomlib_resample_catmull_rom(pts.data(),pts.size(),0,8,
                                        GEOMLIB_CATMULL_CENTRIPETAL,out.data());
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && count==25 &&
                         std::fabs(out[0].x)<1e-4 &&
                         std::fabs(out.back().x-3.0f)<1e-4);
    result.set("count",static_cast<int>(count));
    return result;
}

val meshBilateralDemo() {
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&m,i,i+6,i+5);
        }
    for (size_t i=0;i<m.vertex_count;++i)
        m.vertices[i].z=0.08f*sinf((float)(i*13));
    auto rms=[&]() {
        double s=0.0;
        for (size_t i=0;i<m.vertex_count;++i) s+=m.vertices[i].z*m.vertices[i].z;
        return std::sqrt(s/m.vertex_count);
    };
    double before=rms();
    int rc=geomlib_mesh_bilateral_smooth(&m,nullptr);
    double after=rms();
    val result=val::object();
    result.set("success",rc==GEOMLIB_MESH_SUCCESS && after<before);
    result.set("before",before);
    result.set("after",after);
    geomlib_trimesh_free(&m);
    return result;
}

val meanValueDemo() {
    std::vector<geomlib_vec2> cage = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1)
    };
    std::vector<geomlib_vec2> moved(cage.size());
    for (size_t i = 0; i < cage.size(); ++i)
        moved[i] = geomlib_vec2_create(cage[i].x+2.0f, cage[i].y+1.0f);
    geomlib_vec2 q = geomlib_vec2_create(0.2f,0.3f), out;
    int rc = geomlib_deform_mean_value2d(cage.data(), moved.data(), cage.size(),
                                         &q, 1, &out);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS &&
                          std::fabs(out.x-2.2f) < 1e-3 &&
                          std::fabs(out.y-1.3f) < 1e-3);
    result.set("x", static_cast<double>(out.x));
    result.set("y", static_cast<double>(out.y));
    return result;
}

val abfDemo() {
    geomlib_trimesh m;
    geomlib_trimesh_init(&m, 25, 32);
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 5; ++x)
            geomlib_trimesh_add_vertex(&m,
                geomlib_vec3_create((float)x,(float)y,
                    0.7f*sinf(0.7f*x)*sinf(0.7f*y)), nullptr, nullptr);
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&m,i,i+6,i+5);
        }
    geomlib_uv_coords init = geomlib_uv_coords_create(25,32);
    int rc = geomlib_parameterize_lscm_free(&m,0,4,&init);
    float before = rc == GEOMLIB_MESH_SUCCESS ? geomlib_abf_angle_error(&m,&init) : -1.0f;
    geomlib_abf_params p = geomlib_abf_params_default();
    p.iterations = 30;
    geomlib_uv_coords out = {};
    rc = rc == GEOMLIB_MESH_SUCCESS ? geomlib_parameterize_abf(&m,0,4,&p,&out) : rc;
    float after = rc == GEOMLIB_MESH_SUCCESS ? geomlib_abf_angle_error(&m,&out) : -1.0f;
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && after < before);
    result.set("before", static_cast<double>(before));
    result.set("after", static_cast<double>(after));
    geomlib_uv_coords_free(&out);
    geomlib_uv_coords_free(&init);
    geomlib_trimesh_free(&m);
    return result;
}

val dijkstraDemo() {
    geomlib_trimesh sphere;
    int rc = geomlib_create_sphere(&sphere, 1.0f, 12, 8);
    geomlib_geodesic_distances d = rc == GEOMLIB_MESH_SUCCESS
        ? geomlib_geodesic_create(sphere.vertex_count)
        : geomlib_geodesic_create(0);
    if (rc == GEOMLIB_MESH_SUCCESS)
        rc = geomlib_geodesic_dijkstra_from_vertex(&sphere, 0, &d);
    float mx = 0.0f;
    for (size_t i = 0; i < d.vertex_count; ++i)
        if (d.distances[i] > mx) mx = d.distances[i];
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && mx > 2.0f && mx < 4.0f);
    result.set("max", static_cast<double>(mx));
    geomlib_geodesic_free(&d);
    if (sphere.vertex_count) geomlib_trimesh_free(&sphere);
    return result;
}

val sdfRaycastDemo() {
    geomlib_trimesh sphere;
    geomlib_sdf sdf;
    memset(&sphere, 0, sizeof(sphere));
    memset(&sdf, 0, sizeof(sdf));
    int rc = geomlib_create_sphere(&sphere, 1.0f, 12, 8);
    if (rc == GEOMLIB_MESH_SUCCESS) {
        geomlib_sdf_params p = geomlib_sdf_params_default();
        p.resolution = 24;
        p.compute_sign = true;
        p.padding = 12.0f;
        sdf = geomlib_mesh_to_sdf(&sphere, &p);
    }
    geomlib_sdf_ray_hit hit;
    rc = sdf.values ? geomlib_sdf_raycast(&sdf,
        geomlib_ray_create(geomlib_vec3_create(1.8f,0,0),
                           geomlib_vec3_create(-1,0,0)), 10.0f, 256, &hit)
                    : GEOMLIB_MESH_ERROR_INVALID;
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && hit.hit);
    result.set("t", static_cast<double>(hit.t));
    geomlib_sdf_free(&sdf);
    if (sphere.vertex_count) geomlib_trimesh_free(&sphere);
    return result;
}

val mlsDemo() {
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud, 200);
    for (int i = 0; i < 200; ++i) {
        float x = ((i*37)%101)/100.0f*2.0f-1.0f;
        float y = ((i*73)%97)/100.0f*2.0f-1.0f;
        geomlib_pointcloud_add_point(&cloud, geomlib_vec3_create(x,y,0),
                                     nullptr, 1.0f);
    }
    geomlib_mls_params p = geomlib_mls_params_default();
    p.neighbors = 16;
    geomlib_vec3 out, nrm;
    int rc = geomlib_mls_project(&cloud, geomlib_vec3_create(0.1f,0.2f,0.4f),
                                 &p, &out, &nrm);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && std::fabs(out.z) < 0.03);
    result.set("z", static_cast<double>(out.z));
    geomlib_pointcloud_free(&cloud);
    return result;
}

val simplifyDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i = 0; i < 24; ++i)
        pts.push_back(geomlib_vec3_create((float)i, (i%6==3)?0.5f:0.0f, 0.0f));
    std::vector<uint32_t> idx(pts.size());
    size_t n = 0;
    int rc = geomlib_simplify_polyline(pts.data(), pts.size(), 0.1f,
                                       idx.data(), &n);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && n > 2 && n < pts.size());
    result.set("kept", static_cast<int>(n));
    return result;
}

val polytriHolesDemo() {
    std::vector<geomlib_vec2> outer = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2), geomlib_vec2_create(0,2)
    };
    std::vector<geomlib_vec2> hole = {
        geomlib_vec2_create(0.5f,0.5f), geomlib_vec2_create(1.5f,0.5f),
        geomlib_vec2_create(1.5f,1.5f), geomlib_vec2_create(0.5f,1.5f)
    };
    const geomlib_vec2 *holes[1] = {hole.data()};
    size_t sizes[1] = {hole.size()};
    geomlib_polytri_result r = geomlib_polytri_result_create();
    int rc = geomlib_triangulate_with_holes(outer.data(), outer.size(),
                                            holes, sizes, 1, &r);
    double area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < r.triangle_count; ++i) {
            geomlib_vec2 a=r.vertices[r.triangles[i*3]];
            geomlib_vec2 b=r.vertices[r.triangles[i*3+1]];
            geomlib_vec2 c=r.vertices[r.triangles[i*3+2]];
            area += 0.5*std::fabs((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
        }
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && r.triangle_count > 0);
    result.set("triangles", static_cast<int>(r.triangle_count));
    result.set("area", area);
    geomlib_polytri_result_free(&r);
    return result;
}

val minkowskiDemo() {
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1)
    };
    geomlib_hull2d h = geomlib_hull2d_create();
    int rc = geomlib_minkowski_sum2d(sq.data(), sq.size(), sq.data(), sq.size(), &h);
    double area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < h.count; ++i) {
            size_t j = (i+1)%h.count;
            area += (double)h.vertices[i].x*h.vertices[j].y
                  - (double)h.vertices[j].x*h.vertices[i].y;
        }
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && h.count == 4);
    result.set("count", static_cast<int>(h.count));
    result.set("area", 0.5*std::fabs(area));
    geomlib_hull2d_free(&h);
    return result;
}

val hull2dDemo() {
    std::vector<geomlib_vec2> pts = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2), geomlib_vec2_create(0,2),
        geomlib_vec2_create(1,1)
    };
    geomlib_hull2d h = geomlib_hull2d_create();
    int rc = geomlib_convex_hull2d(pts.data(), pts.size(), &h);
    double area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < h.count; ++i) {
            size_t j = (i+1)%h.count;
            area += (double)h.vertices[i].x*h.vertices[j].y
                  - (double)h.vertices[j].x*h.vertices[i].y;
        }
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && h.count == 4);
    result.set("count", static_cast<int>(h.count));
    result.set("area", 0.5*area);
    geomlib_hull2d_free(&h);
    return result;
}

val nonrigidDemo() {
    geomlib_trimesh src, target;
    int rc = geomlib_create_sphere(&src, 1.0f, 12, 8);
    if (rc == GEOMLIB_MESH_SUCCESS)
        rc = geomlib_trimesh_init(&target, src.vertex_count, src.triangle_count);
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < src.vertex_count; ++i) {
            geomlib_vec3 p = src.vertices[i];
            p.y += 0.3f*sinf(p.x);
            geomlib_trimesh_add_vertex(&target, p, nullptr, nullptr);
        }
        for (size_t i = 0; i < src.triangle_count; ++i)
            geomlib_trimesh_add_triangle(&target, src.indices[i*3],
                src.indices[i*3+1], src.indices[i*3+2]);
    }
    geomlib_nonrigid_result fit = geomlib_nonrigid_result_create();
    if (rc == GEOMLIB_MESH_SUCCESS)
        rc = geomlib_nonrigid_icp(&src, &target, 12, 0.85f, &fit);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS &&
                          fit.final_rms < fit.initial_rms);
    result.set("initial", static_cast<double>(fit.initial_rms));
    result.set("final", static_cast<double>(fit.final_rms));
    geomlib_nonrigid_result_free(&fit);
    if (src.vertex_count) geomlib_trimesh_free(&src);
    if (target.vertex_count) geomlib_trimesh_free(&target);
    return result;
}

val polyboolDemo() {
    std::vector<geomlib_vec2> a = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2), geomlib_vec2_create(0,2)
    };
    std::vector<geomlib_vec2> b = {
        geomlib_vec2_create(1,1), geomlib_vec2_create(3,1),
        geomlib_vec2_create(3,3), geomlib_vec2_create(1,3)
    };
    geomlib_poly_bool_result r = geomlib_poly_bool_result_create();
    int rc = geomlib_polygon_boolean(a.data(), a.size(), b.data(), b.size(),
                                     GEOMLIB_POLY_INTERSECTION, &r);
    double area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS && r.loop_count) {
        for (size_t l = 0; l < r.loop_count; ++l) {
            size_t s=r.offsets[l], e=r.offsets[l+1];
            double v=0.0;
            for (size_t i=s;i<e;++i) {
                size_t j=(i+1<e)?i+1:s;
                v += (double)r.loops[i].x*r.loops[j].y
                   - (double)r.loops[j].x*r.loops[i].y;
            }
            area += std::fabs(0.5*v);
        }
    }
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && std::fabs(area-1.0) < 1e-3);
    result.set("loops", static_cast<int>(r.loop_count));
    result.set("area", area);
    geomlib_poly_bool_result_free(&r);
    return result;
}

val medialAxisDemo() {
    std::vector<geomlib_vec2> poly = {
        geomlib_vec2_create(-1,0), geomlib_vec2_create(0,-1),
        geomlib_vec2_create(1,0), geomlib_vec2_create(0,1)
    };
    geomlib_medial2d m = geomlib_medial2d_create();
    int rc = geomlib_medial_axis2d(poly.data(), poly.size(), 0.08f, &m);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && m.segment_count > 0);
    result.set("segments", static_cast<int>(m.segment_count));
    geomlib_medial2d_free(&m);
    return result;
}

val fpsDemo() {
    geomlib_trimesh sphere;
    int rc = geomlib_create_sphere(&sphere, 1.0f, 12, 8);
    geomlib_fps_result r = geomlib_fps_result_create();
    if (rc == GEOMLIB_MESH_SUCCESS)
        rc = geomlib_farthest_point_sampling(&sphere, 6, 0, &r);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && r.count == 6);
    result.set("count", static_cast<int>(r.count));
    result.set("radius", r.count ? static_cast<double>(r.radii[r.count-1]) : 0.0);
    geomlib_fps_result_free(&r);
    if (rc == GEOMLIB_MESH_SUCCESS || sphere.triangle_count)
        geomlib_trimesh_free(&sphere);
    return result;
}

val offsetDemo() {
    geomlib_trimesh sphere;
    int rc = geomlib_create_sphere(&sphere, 0.5f, 16, 12);
    geomlib_trimesh out;
    memset(&out, 0, sizeof(out));
    geomlib_offset_params p = geomlib_offset_params_default();
    p.resolution = 28;
    if (rc == GEOMLIB_MESH_SUCCESS)
        rc = geomlib_offset_mesh(&sphere, 0.1f, &p, &out);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && out.triangle_count > 0);
    result.set("triangles", static_cast<int>(out.triangle_count));
    result.set("volume", static_cast<double>(
        out.triangle_count ? geomlib_trimesh_volume(&out) : 0.0f));
    geomlib_trimesh_free(&out);
    if (rc == GEOMLIB_MESH_SUCCESS || sphere.triangle_count)
        geomlib_trimesh_free(&sphere);
    return result;
}

val cdt2dDemo() {
    std::vector<geomlib_vec2> pts = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1)
    };
    uint32_t seg[2] = {0,2};
    geomlib_cdt2d c = geomlib_cdt2d_create();
    int rc = geomlib_cdt2d_build(pts.data(), pts.size(), seg, 1, &c);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS &&
                          geomlib_cdt2d_has_edge(&c, 0, 2));
    result.set("triangles", static_cast<int>(c.triangle_count));
    geomlib_cdt2d_free(&c);
    return result;
}

val meshIntersectionDemo() {
    geomlib_trimesh a, b;
    int rc = geomlib_create_box(&a, 1.0f,1.0f,1.0f);
    if (rc == GEOMLIB_MESH_SUCCESS) rc = geomlib_create_box(&b, 1.0f,1.0f,1.0f);
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < b.vertex_count; ++i) {
            b.vertices[i].x += 0.5f;
            b.vertices[i].y += 0.3f;
            b.vertices[i].z += 0.2f;
        }
    geomlib_intersection r = rc == GEOMLIB_MESH_SUCCESS
        ? geomlib_mesh_mesh_intersection(&a, &b) : geomlib_intersection_create();
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && r.segment_count > 0);
    result.set("segments", static_cast<int>(r.segment_count));
    geomlib_intersection_free(&r);
    if (rc == GEOMLIB_MESH_SUCCESS) { geomlib_trimesh_free(&a); geomlib_trimesh_free(&b); }
    return result;
}

val spatialSortDemo() {
    std::vector<geomlib_vec3> pts;
    for (int i = 0; i < 64; ++i)
        pts.push_back(geomlib_vec3_create(((i*37)%101)/100.0f,
                                          ((i*73)%97)/100.0f,
                                          ((i*29)%89)/100.0f));
    std::vector<uint32_t> order(pts.size());
    int rc = geomlib_spatial_sort3(pts.data(), pts.size(), order.data());
    uint64_t code = geomlib_morton3_encode(5, 9, 2);
    uint32_t x=0,y=0,z=0;
    geomlib_morton3_decode(code, &x, &y, &z);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS && x==5 && y==9 && z==2);
    result.set("first", static_cast<int>(order[0]));
    return result;
}

val voxelizeDemo() {
    geomlib_trimesh box;
    int rc = geomlib_create_box(&box, 1.0f, 1.0f, 1.0f);
    geomlib_voxel_grid g = rc == GEOMLIB_MESH_SUCCESS
        ? geomlib_voxelize_mesh(&box, 16, 16, 16) : geomlib_voxel_grid{};
    val result = val::object();
    result.set("success", g.surface != nullptr && g.volume != nullptr);
    result.set("surface", static_cast<int>(g.surface_count));
    result.set("volume", static_cast<int>(g.volume_count));
    geomlib_voxel_grid_free(&g);
    if (rc == GEOMLIB_MESH_SUCCESS) geomlib_trimesh_free(&box);
    return result;
}

val voronoiDemo() {
    std::vector<geomlib_vec2> pts;
    for (int gy = 0; gy < 4; ++gy)
        for (int gx = 0; gx < 4; ++gx)
            pts.push_back(geomlib_vec2_create(0.2f+0.2f*gx, 0.2f+0.2f*gy));
    geomlib_voronoi2d v = geomlib_voronoi2d_create();
    int rc = geomlib_compute_voronoi2d(pts.data(), pts.size(),
        geomlib_vec2_create(0.0f,0.0f), geomlib_vec2_create(1.0f,1.0f), &v);
    double area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < v.site_count; ++i)
            area += geomlib_voronoi2d_cell_area(&v, i);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("sites", static_cast<int>(v.site_count));
    result.set("area", area);
    geomlib_voronoi2d_free(&v);
    return result;
}

val fluidProjectDemo() {
    geomlib_fluid f = geomlib_fluid_create(8,8,8,0.1f);
    for (uint32_t z = 0; z < f.nz; ++z)
        for (uint32_t y = 0; y < f.ny; ++y)
            for (uint32_t x = 0; x < f.nx; ++x) {
                size_t id = x + f.nx*(y + f.ny*z);
                f.u[id] = sinf(0.7f*(float)x);
                f.v[id] = -cosf(0.7f*(float)y);
            }
    auto divRms = [&]() {
        double s = 0.0; int c = 0;
        for (uint32_t z = 1; z+1 < f.nz; ++z)
            for (uint32_t y = 1; y+1 < f.ny; ++y)
                for (uint32_t x = 1; x+1 < f.nx; ++x) {
                    size_t id = x + f.nx*(y + f.ny*z);
                    float d = (f.u[id+1]-f.u[id] + f.v[id+f.nx]-f.v[id])/f.dx;
                    s += d*d; c++;
                }
        return c ? (float)std::sqrt(s/c) : 0.0f;
    };
    float before = divRms();
    int rc = geomlib_fluid_project(&f, 1.0f);
    float after = divRms();
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("before", before);
    result.set("after", after);
    geomlib_fluid_free(&f);
    return result;
}

val hexDominantDemo() {
    geomlib_tetramesh tets;
    geomlib_tetramesh_init(&tets, 8, 5);
    for (int i = 0; i < 8; ++i)
        geomlib_tetramesh_add_vertex(&tets,
            geomlib_vec3_create((float)i, (float)(i%3), (float)(i%2)));
    static const uint32_t t5[5][4] = {
        {0,1,2,3}, {0,1,3,4}, {0,1,4,5}, {0,1,5,6}, {0,1,6,7}
    };
    for (int i = 0; i < 5; ++i)
        geomlib_tetramesh_add_tetrahedron(&tets, t5[i][0], t5[i][1], t5[i][2], t5[i][3]);
    geomlib_hexmesh hm = geomlib_hexmesh_create();
    int rc = geomlib_hex_dominant(&tets, &hm);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("hexes", static_cast<int>(hm.hex_count));
    result.set("tets", static_cast<int>(hm.tet_count));
    geomlib_hexmesh_free(&hm);
    geomlib_tetramesh_free(&tets);
    return result;
}

val sparseDirectSolve(val dense_rows) {
    /* Build a small SPD tridiagonal system from the first n values. */
    size_t len = dense_rows["length"].as<size_t>();
    size_t n = 0;
    while ((n+2)*(n+2) <= len) ++n;
    if (n > 80) n = 80;
    std::vector<uint32_t> rows, cols;
    std::vector<double> vals, b, x(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        rows.push_back((uint32_t)i); cols.push_back((uint32_t)i); vals.push_back(2.0);
        if (i) { rows.push_back((uint32_t)i); cols.push_back((uint32_t)(i-1)); vals.push_back(-1.0); }
        if (i+1<n) { rows.push_back((uint32_t)i); cols.push_back((uint32_t)(i+1)); vals.push_back(-1.0); }
        b.push_back(1.0);
    }
    geomlib_spmat A;
    int rc = n >= 2 ? geomlib_spmat_from_triplets(rows.data(), cols.data(),
        vals.data(), vals.size(), n, n, &A) : GEOMLIB_SPARSE_ERROR_INVALID;
    if (rc == GEOMLIB_SPARSE_SUCCESS) rc = geomlib_solve_sparse_direct(&A, b.data(), x.data(), nullptr);
    val result = val::object();
    result.set("success", rc == GEOMLIB_SPARSE_SUCCESS);
    result.set("size", static_cast<int>(x.size()));
    if (rc == GEOMLIB_SPARSE_SUCCESS) geomlib_spmat_free(&A);
    return result;
}

val octreeQuery(val points, double x, double y, double z, float radius) {
    std::vector<geomlib_vec3> pts;
    size_t n = points["length"].as<size_t>() / 3;
    for (size_t i = 0; i < n; ++i)
        pts.push_back(geomlib_vec3_create(points[i*3].as<float>(),
                                          points[i*3+1].as<float>(),
                                          points[i*3+2].as<float>()));
    geomlib_octree *tree = geomlib_octree_build(pts.data(), pts.size(), nullptr);
    uint32_t idx[256];
    float dst[256];
    size_t count = tree ? geomlib_octree_radius(tree,
        geomlib_vec3_create((float)x,(float)y,(float)z), radius,
        idx, dst, 256) : 0;
    geomlib_octree_free(tree);
    val out = val::object();
    out.set("success", tree != nullptr);
    out.set("count", static_cast<int>(count));
    return out;
}

val bvhRefitDistance(const Mesh& original, const Mesh& moved) {
    geomlib_bvh bvh = {};
    geomlib_bvh_build(&bvh, &original.mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    geomlib_bvh_refit(&bvh, &moved.mesh);
    geomlib_vec3 c;
    float d = 0.0f;
    bool hit = geomlib_bvh_closest_point(&bvh,
        geomlib_vec3_create(0,0,0), &c, &d, NULL);
    geomlib_bvh_free(&bvh);
    val result = val::object();
    result.set("success", hit);
    result.set("distance", d);
    return result;
}

val kdtreeQuery(val points, double x, double y, double z) {
    std::vector<geomlib_vec3> pts;
    size_t n = points["length"].as<size_t>() / 3;
    for (size_t i = 0; i < n; ++i)
        pts.push_back(geomlib_vec3_create(points[i*3].as<float>(),
                                          points[i*3+1].as<float>(),
                                          points[i*3+2].as<float>()));
    geomlib_kdtree *tree = geomlib_kdtree_build(pts.data(), pts.size());
    uint32_t idx[4];
    float dst[4];
    size_t count = tree ? geomlib_kdtree_knn(tree,
        geomlib_vec3_create((float)x,(float)y,(float)z), 4, idx, dst) : 0;
    geomlib_kdtree_free(tree);
    val out = val::object();
    out.set("success", count == 4);
    out.set("count", static_cast<int>(count));
    out.set("distance", count ? dst[0] : 0.0f);
    return out;
}

val shapeDiameter(const Mesh& mesh) {
    geomlib_shape_diameter d = geomlib_shape_diameter_create(0);
    int rc = geomlib_compute_shape_diameter(&mesh.mesh, 8, &d);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("count", static_cast<int>(d.count));
    val values = val::array();
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < d.count; ++i) values.call<void>("push", d.values[i]);
    }
    result.set("values", values);
    geomlib_shape_diameter_free(&d);
    return result;
}

val crossField(const Mesh& mesh) {
    geomlib_cross_field f = geomlib_cross_field_create(0);
    int rc = geomlib_compute_cross_field(&mesh.mesh, &f);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("count", static_cast<int>(f.count));
    val angles = val::array();
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < f.count; ++i) angles.call<void>("push", f.angles[i]);
    }
    result.set("angles", angles);
    geomlib_cross_field_free(&f);
    return result;
}

val skinRotateDualQuat(const Mesh& mesh, double angle) {
    Mesh* out = new Mesh();
    std::vector<geomlib_skin_weight> weights(mesh.mesh.vertex_count);
    for (auto &w : weights) { w.bone_indices[0]=0; w.weights[0]=1.0f; }
    geomlib_mat4 rot = geomlib_mat4_identity();
    rot.m[0] = (float)cos(angle); rot.m[1] = (float)sin(angle);
    rot.m[4] = (float)-sin(angle); rot.m[5] = (float)cos(angle);
    int rc = geomlib_skin_mesh_dualquat(&mesh.mesh, weights.data(), &rot, 1, &out->mesh);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("vertices", static_cast<int>(out->mesh.vertex_count));
    delete out;
    return result;
}

val skinRotate(const Mesh& mesh, double angle) {
    Mesh* out = new Mesh();
    std::vector<geomlib_skin_weight> weights(mesh.mesh.vertex_count);
    for (auto &w : weights) { w.bone_indices[0]=0; w.weights[0]=1.0f; }
    geomlib_mat4 rot = geomlib_mat4_identity();
    rot.m[0] = (float)cos(angle); rot.m[1] = (float)sin(angle);
    rot.m[4] = (float)-sin(angle); rot.m[5] = (float)cos(angle);
    int rc = geomlib_skin_mesh(&mesh.mesh, weights.data(), &rot, 1, &out->mesh);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("vertices", static_cast<int>(out->mesh.vertex_count));
    delete out;
    return result;
}

val geodesicUnfold(const Mesh& mesh, int source, int target) {
    double length = 0.0;
    int rc = geomlib_geodesic_length_unfolding(&mesh.mesh, (uint32_t)source,
                                               (uint32_t)target, &length);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("length", length);
    return result;
}

val geodesicFMM(const Mesh& mesh, int source) {
    geomlib_geodesic_fmm_result r = geomlib_geodesic_fmm_result_create(0);
    int rc = geomlib_geodesic_fast_marching(&mesh.mesh, (uint32_t)source, 4, &r);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("count", static_cast<int>(r.count));
    val values = val::array();
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < r.count; ++i) values.call<void>("push", r.distances[i]);
    }
    result.set("distances", values);
    geomlib_geodesic_fmm_result_free(&r);
    return result;
}

val delaunay2d(val points) {
    size_t n = points["length"].as<size_t>() / 2;
    std::vector<geomlib_vec2> pts;
    for (size_t i = 0; i < n; ++i)
        pts.push_back(geomlib_vec2_create(points[i*2].as<float>(),
                                          points[i*2+1].as<float>()));
    geomlib_delaunay2d_result r = geomlib_delaunay2d_result_create();
    int rc = n >= 3 ? geomlib_delaunay2d(pts.data(), n, &r) : GEOMLIB_MESH_ERROR_INVALID;
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("triangles", static_cast<int>(r.triangle_count));
    geomlib_delaunay2d_result_free(&r);
    return result;
}

val marchingSquares(val values, float h) {
    size_t n = values["length"].as<size_t>();
    geomlib_contour_result r = geomlib_contour_result_create();
    std::vector<float> grid;
    for (size_t i = 0; i < n; ++i) grid.push_back(values[i].as<float>());
    size_t side = 0;
    while ((side+1)*(side+1) <= n) ++side;
    int rc = side*side == n && side >= 2
        ? geomlib_marching_squares(grid.data(), side, side, h,
                                   geomlib_vec2_create(0,0), &r)
        : GEOMLIB_MESH_ERROR_INVALID;
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("contours", static_cast<int>(r.polyline_count));
    result.set("points", static_cast<int>(r.point_count));
    geomlib_contour_result_free(&r);
    return result;
}

val clusterDecimate(const Mesh& mesh, float voxel_size) {
    Mesh* out = new Mesh();
    int rc = geomlib_decimate_vertex_clustering(&mesh.mesh, voxel_size, &out->mesh);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("vertices", static_cast<int>(out->mesh.vertex_count));
    result.set("triangles", static_cast<int>(out->mesh.triangle_count));
    delete out;
    return result;
}

val partitionMesh(const Mesh& mesh, int parts) {
    geomlib_partition_result r = geomlib_partition_result_create();
    int rc = geomlib_partition_mesh(&mesh.mesh, (size_t)parts, &r);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    val sizes = val::array();
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < r.part_count; ++i)
            sizes.call<void>("push", static_cast<int>(r.part_sizes[i]));
    }
    result.set("sizes", sizes);
    geomlib_partition_result_free(&r);
    return result;
}

val binaryPlyRoundtrip(const Mesh& mesh) {
    void *buf = nullptr;
    size_t size = 0;
    int rc = geomlib_mesh_save_mem(&mesh.mesh, "plyb", &buf, &size);
    Mesh* loaded = new Mesh();
    if (rc == GEOMLIB_IO_SUCCESS) {
        rc = geomlib_mesh_load_mem(&loaded->mesh, "plyb", buf, size);
        free(buf);
    }
    val result = val::object();
    result.set("success", rc == GEOMLIB_IO_SUCCESS);
    result.set("vertices", static_cast<int>(loaded->mesh.vertex_count));
    result.set("triangles", static_cast<int>(loaded->mesh.triangle_count));
    delete loaded;
    return result;
}

val fastWinding(const Mesh& mesh, double x, double y, double z) {
    /* The tree borrows the mesh, so keep it alive for the query. */
    geomlib_winding_tree *tree = geomlib_winding_tree_build(&mesh.mesh, 24);
    float value = tree ? geomlib_fast_winding_number(
        tree, geomlib_vec3_create((float)x, (float)y, (float)z)) : 0.0f;
    geomlib_winding_tree_free(tree);
    val result = val::object();
    result.set("value", value);
    return result;
}

val fairMesh(const Mesh& mesh, float lambda, int iterations) {
    Mesh* out = new Mesh();
    geomlib_trimesh_init(&out->mesh, mesh.mesh.vertex_count, mesh.mesh.triangle_count);
    memcpy(out->mesh.vertices, mesh.mesh.vertices,
           mesh.mesh.vertex_count * sizeof(geomlib_vec3));
    memcpy(out->mesh.indices, mesh.mesh.indices,
           mesh.mesh.triangle_count * 3 * sizeof(uint32_t));
    out->mesh.vertex_count = mesh.mesh.vertex_count;
    out->mesh.triangle_count = mesh.mesh.triangle_count;
    geomlib_fairing_params params = geomlib_fairing_params_default();
    params.lambda = lambda;
    params.iterations = iterations;
    int rc = geomlib_fair_mesh(&out->mesh, &params);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("vertices", static_cast<int>(out->mesh.vertex_count));
    delete out;
    return result;
}

val convexDecompose(const Mesh& mesh, int pieces) {
    geomlib_convex_decomposition d = geomlib_convex_decomposition_create();
    int rc = geomlib_approx_convex_decomposition(&mesh.mesh, (size_t)pieces, &d);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("pieces", static_cast<int>(d.piece_count));
    int triangles = 0;
    for (size_t i = 0; i < d.piece_count; ++i)
        triangles += (int)d.pieces[i].triangle_count;
    result.set("triangles", triangles);
    geomlib_convex_decomposition_free(&d);
    return result;
}

val sliceMesh(const Mesh& mesh, double nx, double ny, double nz, double d) {
    geomlib_slice_result r = geomlib_slice_result_create();
    geomlib_vec3 n = geomlib_vec3_create((float)nx, (float)ny, (float)nz);
    int rc = geomlib_slice_mesh(&mesh.mesh, geomlib_vec3_mul(n, (float)d), n,
                                1e-4f, &r);
    val out = val::object();
    out.set("success", rc == GEOMLIB_MESH_SUCCESS);
    out.set("polylines", static_cast<int>(r.polyline_count));
    out.set("points", static_cast<int>(r.point_count));
    geomlib_slice_result_free(&r);
    return out;
}

Mesh* alphaWrap(const Mesh& mesh, float alpha, int resolution) {
    geomlib_alpha_wrap_params params = geomlib_alpha_wrap_params_default();
    params.alpha = alpha;
    params.resolution = (uint32_t)resolution;
    Mesh* out = new Mesh();
    if (geomlib_alpha_wrap(&mesh.mesh, &params, &out->mesh) != GEOMLIB_MESH_SUCCESS) {
        delete out;
        return nullptr;
    }
    return out;
}

val denoisePointCloud(const PointCloud& cloud, int iterations) {
    PointCloud* copy = new PointCloud();
    for (size_t i = 0; i < cloud.pc.point_count; ++i) {
        geomlib_vec3 n = cloud.pc.normals ? cloud.pc.normals[i] : geomlib_vec3_create(0,0,1);
        copy->addPointWithNormal(cloud.pc.points[i].x, cloud.pc.points[i].y,
                                 cloud.pc.points[i].z, n.x, n.y, n.z);
    }
    geomlib_bilateral_params params = geomlib_bilateral_params_default();
    params.iterations = iterations;
    int rc = geomlib_pointcloud_bilateral_smooth(&copy->pc, &params);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("count", static_cast<int>(copy->pc.point_count));
    delete copy;
    return result;
}

val vsaSimplify(const Mesh& mesh, int regions) {
    geomlib_vsa_result result = geomlib_vsa_result_create();
    Mesh* out = new Mesh();
    int rc = geomlib_vsa_simplify(&mesh.mesh, (size_t)regions, 5, &result, &out->mesh);
    val v = val::object();
    v.set("success", rc == GEOMLIB_MESH_SUCCESS);
    v.set("regions", static_cast<int>(result.region_count));
    v.set("triangles", static_cast<int>(out->mesh.triangle_count));
    geomlib_vsa_result_free(&result);
    delete out;
    return v;
}

val offsetPolygon(val points, float distance) {
    geomlib_polygon2d poly = {};
    size_t n = points["length"].as<size_t>() / 2;
    if (n >= 3) {
        geomlib_polygon2d_init(&poly, n);
        for (size_t i = 0; i < n; ++i) {
            geomlib_vec2 p = { points[i*2].as<float>(), points[i*2+1].as<float>() };
            geomlib_polygon2d_add_vertex(&poly, p);
        }
    }
    geomlib_offset_result result = geomlib_offset_result_create();
    int rc = n >= 3 ? geomlib_offset_polygon(&poly, distance, 0, &result)
                    : GEOMLIB_MESH_ERROR_INVALID;
    val out = val::object();
    out.set("success", rc == GEOMLIB_MESH_SUCCESS);
    out.set("loops", static_cast<int>(result.count));
    out.set("vertices", static_cast<int>(result.count ? result.polygons[0].vertex_count : 0));
    geomlib_offset_result_free(&result);
    geomlib_polygon2d_free(&poly);
    return out;
}

val straightSkeleton(val points) {
    geomlib_polygon2d poly = {};
    size_t n = points["length"].as<size_t>() / 2;
    if (n >= 3) {
        geomlib_polygon2d_init(&poly, n);
        for (size_t i = 0; i < n; ++i) {
            geomlib_vec2 p = {
                points[i*2].as<float>(),
                points[i*2+1].as<float>()
            };
            geomlib_polygon2d_add_vertex(&poly, p);
        }
    }
    geomlib_straight_skeleton s = geomlib_straight_skeleton_create();
    int rc = n >= 3 ? geomlib_compute_straight_skeleton(&poly, &s) : GEOMLIB_MESH_ERROR_INVALID;
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("nodes", static_cast<int>(s.point_count));
    result.set("edges", static_cast<int>(s.edge_count));
    geomlib_straight_skeleton_free(&s);
    geomlib_polygon2d_free(&poly);
    return result;
}

val extractSkeleton(const Mesh& mesh, float ratio) {
    geomlib_skeleton_params params = geomlib_skeleton_params_default();
    params.collapse_ratio = ratio;
    params.max_iterations = 60;
    geomlib_curve_skeleton s = geomlib_curve_skeleton_create();
    int rc = geomlib_extract_curve_skeleton(&mesh.mesh, &params, &s);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    result.set("nodes", static_cast<int>(s.point_count));
    result.set("edges", static_cast<int>(s.edge_count));
    geomlib_curve_skeleton_free(&s);
    return result;
}

val computeSaliency(const Mesh& mesh) {
    geomlib_saliency_result r = geomlib_saliency_result_create(mesh.mesh.vertex_count);
    int rc = geomlib_compute_saliency(&mesh.mesh, 3, &r);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    val values = val::array();
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < r.vertex_count; ++i)
            values.call<void>("push", r.values[i]);
    }
    result.set("values", values);
    geomlib_saliency_result_free(&r);
    return result;
}

Mesh* dualContouring(const Mesh& sdf_mesh, int resolution) {
    geomlib_sdf sdf = geomlib_mesh_to_sdf_resolution(&sdf_mesh.mesh, (uint32_t)resolution);
    Mesh* out = new Mesh();
    if (!sdf.values || geomlib_sdf_dual_contouring(&sdf, &out->mesh) != GEOMLIB_MESH_SUCCESS) {
        delete out;
        out = nullptr;
    }
    geomlib_sdf_free(&sdf);
    return out;
}

Mesh* marchingTets(const Mesh& sdf_mesh, int resolution) {
    geomlib_sdf sdf = geomlib_mesh_to_sdf_resolution(&sdf_mesh.mesh, (uint32_t)resolution);
    Mesh* out = new Mesh();
    if (!sdf.values || geomlib_sdf_marching_tets(&sdf, &out->mesh) != GEOMLIB_MESH_SUCCESS) {
        delete out;
        out = nullptr;
    }
    geomlib_sdf_free(&sdf);
    return out;
}

val remeshAnisotropic(const Mesh& mesh, float target_length, int iterations) {
    Mesh* out = new Mesh();
    geomlib_trimesh_init(&out->mesh, mesh.mesh.vertex_count, mesh.mesh.triangle_count);
    memcpy(out->mesh.vertices, mesh.mesh.vertices,
           mesh.mesh.vertex_count * sizeof(geomlib_vec3));
    memcpy(out->mesh.indices, mesh.mesh.indices,
           mesh.mesh.triangle_count * 3 * sizeof(uint32_t));
    out->mesh.vertex_count = mesh.mesh.vertex_count;
    out->mesh.triangle_count = mesh.mesh.triangle_count;
    geomlib_anisotropic_params params = geomlib_anisotropic_params_default(target_length);
    params.iterations = iterations;
    int rc = geomlib_remesh_anisotropic(&out->mesh, &params);
    if (rc != GEOMLIB_MESH_SUCCESS) { delete out; return val::undefined(); }
    val result = val::object();
    result.set("vertices", static_cast<int>(out->mesh.vertex_count));
    result.set("triangles", static_cast<int>(out->mesh.triangle_count));
    delete out;
    return result;
}

val remeshToQuads(const Mesh& mesh) {
    geomlib_quadmesh result = {};
    int rc = geomlib_remesh_to_quads(&mesh.mesh, nullptr, &result);
    val out = val::object();
    out.set("success", rc == GEOMLIB_MESH_SUCCESS);
    out.set("quads", static_cast<int>(result.quad_count));
    out.set("triangles", static_cast<int>(result.triangle_count));
    geomlib_quadmesh_free(&result);
    return out;
}

val detectPlanes(const PointCloud& cloud, float threshold, int max_shapes) {
    size_t count = 0;
    geomlib_detected_shape *shapes = geomlib_detect_shapes_ransac(
        &cloud.pc, GEOMLIB_SHAPE_PLANE, threshold, (size_t)max_shapes, 42, &count);
    val result = val::object();
    result.set("count", static_cast<int>(count));
    val arr = val::array();
    for (size_t i = 0; i < count; ++i) {
        val s = val::object();
        s.set("inliers", static_cast<int>(shapes[i].inlier_count));
        s.set("nx", shapes[i].normal.x);
        s.set("ny", shapes[i].normal.y);
        s.set("nz", shapes[i].normal.z);
        arr.call<void>("push", s);
    }
    result.set("shapes", arr);
    geomlib_detected_shapes_free(shapes, count);
    return result;
}

val manifoldHarmonics(const Mesh& mesh, int count) {
    geomlib_harmonics_result r = {};
    int rc = geomlib_compute_manifold_harmonics(&mesh.mesh, (size_t)count, &r);
    val result = val::object();
    result.set("success", rc == GEOMLIB_MESH_SUCCESS);
    val eigen = val::array();
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < r.count; ++i) eigen.call<void>("push", r.eigenvalues[i]);
    }
    result.set("eigenvalues", eigen);
    geomlib_harmonics_result_free(&r);
    return result;
}

Mesh* deformARAP(const Mesh& mesh, int handle_vertex, double tx, double ty, double tz) {
    geomlib_deform_handle handle;
    handle.vertex = (uint32_t)handle_vertex;
    handle.target = geomlib_vec3_create((float)tx, (float)ty, (float)tz);
    Mesh* out = new Mesh();
    if (geomlib_deform_arap(&mesh.mesh, &handle, 1, 5, &out->mesh) != GEOMLIB_MESH_SUCCESS) {
        delete out;
        return nullptr;
    }
    return out;
}

Mesh* remeshACVD(const Mesh& mesh, int target_vertices) {
    Mesh* out = new Mesh();
    geomlib_acvd_params params = geomlib_acvd_params_default((uint32_t)target_vertices);
    params.lloyd_iterations = 3;
    int rc = geomlib_remesh_acvd(&mesh.mesh, &params, &out->mesh);
    if (rc != GEOMLIB_MESH_SUCCESS) { delete out; return nullptr; }
    return out;
}

Mesh* reconstructPoisson(const PointCloud& cloud) {
    Mesh* out = new Mesh();
    int rc = geomlib_reconstruct_poisson(&cloud.pc, nullptr, &out->mesh);
    if (rc != GEOMLIB_MESH_SUCCESS) { delete out; return nullptr; }
    return out;
}

val blueNoiseSurface(const Mesh& mesh, float min_distance) {
    geomlib_sample_set samples = {};
    size_t n = geomlib_blue_noise_surface(&mesh.mesh, min_distance, nullptr, &samples);
    val result = val::object();
    result.set("count", static_cast<int>(n));
    val pts = val::array();
    for (size_t i = 0; i < samples.count; i++) {
        pts.call<void>("push", samples.points[i].x);
        pts.call<void>("push", samples.points[i].y);
        pts.call<void>("push", samples.points[i].z);
    }
    result.set("points", pts);
    geomlib_sample_set_free(&samples);
    return result;
}

val rigidICP(const Mesh& source, const Mesh& target, int max_iterations) {
    geomlib_icp_params params = geomlib_icp_params_default();
    params.max_iterations = max_iterations;
    geomlib_mat4 transform = geomlib_mat4_identity();
    float rms = geomlib_icp_rigid(&source.mesh, &target.mesh, &params, &transform);
    val result = val::object();
    result.set("rms", rms);
    val m = val::array();
    for (int i = 0; i < 16; ++i) m.call<void>("push", transform.m[i]);
    result.set("transform", m);
    return result;
}

Mesh* meshBoolean(const Mesh& a, const Mesh& b, int op) {
    Mesh* out = new Mesh();
    int rc = geomlib_mesh_boolean(&a.mesh, &b.mesh,
                                  static_cast<geomlib_boolean_op>(op), 0.0f, &out->mesh);
    if (rc != GEOMLIB_MESH_SUCCESS) { delete out; return nullptr; }
    return out;
}

Mesh* subdivideCatmullClark(const Mesh& mesh) {
    Mesh* out = new Mesh();
    if (geomlib_subdivide_catmull_clark(&mesh.mesh, &out->mesh) != GEOMLIB_MESH_SUCCESS) {
        delete out;
        return nullptr;
    }
    return out;
}

// ============================================================================
// Memory-buffer I/O and robust predicates (new foundational modules)
// ============================================================================

Mesh* loadMeshFromString(const std::string& text, const std::string& format) {
    Mesh* m = new Mesh();
    geomlib_trimesh_free(&m->mesh); // load_mem initializes the mesh itself
    const char* fmt = format.empty() ? nullptr : format.c_str();
    if (geomlib_mesh_load_mem(&m->mesh, fmt, text.data(), text.size()) != GEOMLIB_IO_SUCCESS) {
        delete m;
        return nullptr;
    }
    return m;
}

std::string saveMeshToString(const Mesh& mesh, const std::string& format) {
    void* buf = nullptr;
    size_t size = 0;
    if (geomlib_mesh_save_mem(&mesh.mesh, format.c_str(), &buf, &size) != GEOMLIB_IO_SUCCESS) {
        return std::string();
    }
    std::string result((const char*)buf, size);
    free(buf);
    return result;
}

int predOrient2d(double ax, double ay, double bx, double by, double cx, double cy) {
    return geomlib_pred_orient2d(ax, ay, bx, by, cx, cy);
}

int predOrient3d(double ax, double ay, double az, double bx, double by, double bz,
                 double cx, double cy, double cz, double dx, double dy, double dz) {
    return geomlib_pred_orient3d(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
}

int predInCircle(double ax, double ay, double bx, double by,
                 double cx, double cy, double dx, double dy) {
    return geomlib_pred_incircle(ax, ay, bx, by, cx, cy, dx, dy);
}

int predInSphere(double ax, double ay, double az, double bx, double by, double bz,
                 double cx, double cy, double cz, double dx, double dy, double dz,
                 double ex, double ey, double ez) {
    return geomlib_pred_insphere(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz,
                                 ex, ey, ez);
}

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
        .function("getUVs", &Mesh::getUVs)
        ;

    // BVH
    class_<BVH>("BVH")
        .constructor<const Mesh*>()
        .function("traceRay", &BVH::traceRay)
        .function("getStats", &BVH::getStats)
        .function("getBoundingBoxes", &BVH::getBoundingBoxes, allow_raw_pointers())
        ;

    // Wave 1 parameterization and Catmull-Clark
    function("parameterizeMesh", &parameterizeMesh, allow_raw_pointers());
    function("subdivideCatmullClark", &subdivideCatmullClark, allow_raw_pointers());
    function("meshBoolean", &meshBoolean, allow_raw_pointers());
    function("rigidICP", &rigidICP, allow_raw_pointers());
    function("blueNoiseSurface", &blueNoiseSurface, allow_raw_pointers());
    function("reconstructPoisson", &reconstructPoisson, allow_raw_pointers());
    function("remeshACVD", &remeshACVD, allow_raw_pointers());
    function("deformARAP", &deformARAP, allow_raw_pointers());
    function("manifoldHarmonics", &manifoldHarmonics, allow_raw_pointers());
    function("detectPlanes", &detectPlanes, allow_raw_pointers());
    function("remeshToQuads", &remeshToQuads, allow_raw_pointers());
    function("remeshAnisotropic", &remeshAnisotropic, allow_raw_pointers());
    function("marchingTets", &marchingTets, allow_raw_pointers());
    function("dualContouring", &dualContouring, allow_raw_pointers());
    function("computeSaliency", &computeSaliency, allow_raw_pointers());
    function("extractSkeleton", &extractSkeleton, allow_raw_pointers());
    function("straightSkeleton", &straightSkeleton);
    function("offsetPolygon", &offsetPolygon);
    function("vsaSimplify", &vsaSimplify, allow_raw_pointers());
    function("denoisePointCloud", &denoisePointCloud, allow_raw_pointers());
    function("alphaWrap", &alphaWrap, allow_raw_pointers());
    function("sliceMesh", &sliceMesh, allow_raw_pointers());
    function("convexDecompose", &convexDecompose, allow_raw_pointers());
    function("fairMesh", &fairMesh, allow_raw_pointers());
    function("fastWinding", &fastWinding, allow_raw_pointers());
    function("binaryPlyRoundtrip", &binaryPlyRoundtrip, allow_raw_pointers());
    function("partitionMesh", &partitionMesh, allow_raw_pointers());
    function("clusterDecimate", &clusterDecimate, allow_raw_pointers());
    function("marchingSquares", &marchingSquares);
    function("delaunay2d", &delaunay2d);
    function("geodesicFMM", &geodesicFMM, allow_raw_pointers());
    function("geodesicUnfold", &geodesicUnfold, allow_raw_pointers());
    function("skinRotate", &skinRotate, allow_raw_pointers());
    function("skinRotateDualQuat", &skinRotateDualQuat, allow_raw_pointers());
    function("crossField", &crossField, allow_raw_pointers());
    function("shapeDiameter", &shapeDiameter, allow_raw_pointers());
    function("kdtreeQuery", &kdtreeQuery);
    function("bvhRefitDistance", &bvhRefitDistance, allow_raw_pointers());
    function("octreeQuery", &octreeQuery);
    function("sparseDirectSolve", &sparseDirectSolve);
    function("hexDominantDemo", &hexDominantDemo);
    function("fluidProjectDemo", &fluidProjectDemo);
    function("voronoiDemo", &voronoiDemo);
    function("voxelizeDemo", &voxelizeDemo);
    function("spatialSortDemo", &spatialSortDemo);
    function("meshIntersectionDemo", &meshIntersectionDemo);
    function("cdt2dDemo", &cdt2dDemo);
    function("offsetDemo", &offsetDemo);
    function("fpsDemo", &fpsDemo);
    function("medialAxisDemo", &medialAxisDemo);
    function("polyboolDemo", &polyboolDemo);
    function("nonrigidDemo", &nonrigidDemo);
    function("hull2dDemo", &hull2dDemo);
    function("minkowskiDemo", &minkowskiDemo);
    function("polytriHolesDemo", &polytriHolesDemo);
    function("simplifyDemo", &simplifyDemo);
    function("mlsDemo", &mlsDemo);
    function("sdfRaycastDemo", &sdfRaycastDemo);
    function("dijkstraDemo", &dijkstraDemo);
    function("abfDemo", &abfDemo);
    function("meanValueDemo", &meanValueDemo);
    function("meshBilateralDemo", &meshBilateralDemo);
    function("splineDemo", &splineDemo);
    function("qualityFlipDemo", &qualityFlipDemo);
    function("refineDemo", &refineDemo);
    function("collapseDemo", &collapseDemo);
    function("refine2dDemo", &refine2dDemo);
    function("surfaceVoronoiDemo", &surfaceVoronoiDemo);
    function("icpPlaneDemo", &icpPlaneDemo);
    function("meshContourDemo", &meshContourDemo);
    function("regionGrowingDemo", &regionGrowingDemo);
    function("pointCurvatureDemo", &pointCurvatureDemo);
    function("hksDemo", &hksDemo);
    function("convexDecomposeDemo", &convexDecomposeDemo);
    function("boundingSphereDemo", &boundingSphereDemo);
    function("miniballDemo", &miniballDemo);
    function("massPropertiesDemo", &massPropertiesDemo);
    function("inscribedCircleDemo", &inscribedCircleDemo);
    function("minAreaRectDemo", &minAreaRectDemo);
    function("polylineIntersectionDemo", &polylineIntersectionDemo);
    function("ambientOcclusionDemo", &ambientOcclusionDemo);
    function("spectralSegmentDemo", &spectralSegmentDemo);
    function("meshDistanceDemo", &meshDistanceDemo);
    function("symmetryDemo", &symmetryDemo);
    function("kmeansDemo", &kmeansDemo);
    function("dbscanDemo", &dbscanDemo);
    function("resampleDemo", &resampleDemo);
    function("polygonMomentsDemo", &polygonMomentsDemo);
    function("polygonClipDemo", &polygonClipDemo);
    function("minCircleDemo", &minCircleDemo);
    function("geodesicPathDemo", &geodesicPathDemo);
    function("geodesicBallDemo", &geodesicBallDemo);
    function("rigidFitDemo", &rigidFitDemo);
    function("circleFitDemo", &circleFitDemo);
    function("sphereFitDemo", &sphereFitDemo);
    function("planeFitDemo", &planeFitDemo);
    function("chaikinDemo", &chaikinDemo);
    function("ellipseFitDemo", &ellipseFitDemo);
    function("cylinderFitDemo", &cylinderFitDemo);
    function("segmentDistanceDemo", &segmentDistanceDemo);
    function("rayPlaneDemo", &rayPlaneDemo);
    function("triTriDemo", &triTriDemo);
    function("polyboolHolesDemo", &polyboolHolesDemo);
    function("geodesicDiameterDemo", &geodesicDiameterDemo);
    function("cage3dDemo", &cage3dDemo);
    function("poisson3dDemo", &poisson3dDemo);
    function("pointInPolygonDemo", &pointInPolygonDemo);
    function("laplacianDemo", &laplacianDemo);
    function("quaternionDemo", &quaternionDemo);
    function("meshClipDemo", &meshClipDemo);
    function("segment2dDemo", &segment2dDemo);
    function("mstDemo", &mstDemo);
    function("pcaDemo", &pcaDemo);
    function("sphereIntersectionDemo", &sphereIntersectionDemo);
    function("frechetDemo", &frechetDemo);
    function("poisson2dDemo", &poisson2dDemo);
    function("curveSmoothDemo", &curveSmoothDemo);
    function("cage2dHarmonicDemo", &cage2dHarmonicDemo);
    function("raster2dDemo", &raster2dDemo);
    function("visibilityDemo", &visibilityDemo);
    function("meshOverlapDemo", &meshOverlapDemo);
    function("solidAngleDemo", &solidAngleDemo);
    function("vwDemo", &vwDemo);
    function("bezierDemo", &bezierDemo);
    function("raySphereDemo", &raySphereDemo);
    function("planePlaneDemo", &planePlaneDemo);
    function("sphereAabbDemo", &sphereAabbDemo);
    function("circleCircleDemo", &circleCircleDemo);
    function("segmentTriangleDemo", &segmentTriangleDemo);
    function("rayCylinderDemo", &rayCylinderDemo);
    function("capsuleDemo", &capsuleDemo);
    function("sphereTriangleDemo", &sphereTriangleDemo);
    function("triAabbDemo", &triAabbDemo);
    function("polygonDiameterDemo", &polygonDiameterDemo);
    function("polygonWidthDemo", &polygonWidthDemo);
    function("aabbDistanceDemo", &aabbDistanceDemo);
    function("rayObbDemo", &rayObbDemo);
    function("triTriDistanceDemo", &triTriDistanceDemo);
    function("planeSphereDemo", &planeSphereDemo);
    function("circleLineDemo", &circleLineDemo);
    function("rayTetraDemo", &rayTetraDemo);
    function("polylineClosestDemo", &polylineClosestDemo);
    function("capsuleAabbDemo", &capsuleAabbDemo);
    function("polygonDistance2dDemo", &polygonDistance2dDemo);
    function("triTetraDemo", &triTetraDemo);
    function("sphereCapsuleDemo", &sphereCapsuleDemo);
    function("tetraTetraDemo", &tetraTetraDemo);
    function("sphereObbDemo", &sphereObbDemo);
    function("circleAabbDemo", &circleAabbDemo);
    function("capsuleObbDemo", &capsuleObbDemo);
    function("rayConeDemo", &rayConeDemo);
    function("triObbDemo", &triObbDemo);
    function("aabbObbDemo", &aabbObbDemo);
    function("obbObbDemo", &obbObbDemo);
    function("ellipseDemo", &ellipseDemo);
    function("capsuleTriangleDemo", &capsuleTriangleDemo);
    function("ellipseClosestDemo", &ellipseClosestDemo);
    function("ellipseMeasuresDemo", &ellipseMeasuresDemo);
    function("circleTriangleDemo", &circleTriangleDemo);
    function("line2dDemo", &line2dDemo);
    function("arcLengthDemo", &arcLengthDemo);
    function("circleAreaDemo", &circleAreaDemo);
    function("polygonIntersectionAreaDemo", &polygonIntersectionAreaDemo);
    function("circleEllipseDemo", &circleEllipseDemo);
    function("geodesicMmpDemo", &geodesicMmpDemo);
    function("metricDemo", &metricDemo);
    function("reconstructRBF", &reconstructRBF);

    // Memory-buffer I/O and robust predicates
    function("loadMeshFromString", &loadMeshFromString, allow_raw_pointers());
    function("saveMeshToString", &saveMeshToString);
    function("predOrient2d", &predOrient2d);
    function("predOrient3d", &predOrient3d);
    function("predInCircle", &predInCircle);
    function("predInSphere", &predInSphere);

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
        .function("extractSurface", &SDF::extractSurface, return_value_policy::take_ownership())
        .function("getStats", &SDF::getStats)
        ;

    function("generateSDF", &generateSDF, allow_raw_pointers());
    function("generateSDFResolution", &generateSDFResolution, allow_raw_pointers());

    // Subdivision
    function("subdivideLoop", &subdivideLoop, allow_raw_pointers());
    function("subdivideButterfly", &subdivideButterfly, allow_raw_pointers());
    function("subdivideMidpoint", &subdivideMidpoint, allow_raw_pointers());
}
