/*
 * geomlib - Monte Carlo Geometry Processing Tests
 * Tests for Monte Carlo sampling and integration
 */

#include <geomlib/geomlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define TEST_EPSILON 1e-5f
#define TEST_EPSILON_LOOSE 1e-3f

void test_triangle_sampling_uniform() {
    printf("Testing uniform triangle sampling...\n");

    geomlib_vec3 v0 = {0.0f, 0.0f, 0.0f};
    geomlib_vec3 v1 = {1.0f, 0.0f, 0.0f};
    geomlib_vec3 v2 = {0.0f, 1.0f, 0.0f};

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 10000;
    geomlib_vec3 centroid = {0, 0, 0};

    for (int i = 0; i < num_samples; i++) {
        geomlib_vec3 sample = geomlib_sample_triangle_uniform(v0, v1, v2, &rng);

        /* Check sample is on the plane z=0 */
        assert(fabsf(sample.z) < TEST_EPSILON);

        /* Check sample is inside triangle (barycentric coordinates >= 0) */
        assert(sample.x >= -TEST_EPSILON);
        assert(sample.y >= -TEST_EPSILON);
        assert(sample.x + sample.y <= 1.0f + TEST_EPSILON);

        centroid.x += sample.x;
        centroid.y += sample.y;
        centroid.z += sample.z;
    }

    centroid.x /= num_samples;
    centroid.y /= num_samples;
    centroid.z /= num_samples;

    /* Expected centroid of triangle (0,0,0), (1,0,0), (0,1,0) is (1/3, 1/3, 0) */
    float expected_x = 1.0f / 3.0f;
    float expected_y = 1.0f / 3.0f;

    assert(fabsf(centroid.x - expected_x) < 0.01f);
    assert(fabsf(centroid.y - expected_y) < 0.01f);
    assert(fabsf(centroid.z) < 0.01f);

    printf("  ✓ Uniform sampling produces correct distribution (centroid: %.4f, %.4f, %.4f)\n",
           centroid.x, centroid.y, centroid.z);
}

void test_triangle_sampling_barycentric() {
    printf("Testing triangle sampling with barycentric coordinates...\n");

    geomlib_vec3 v0 = {0.0f, 0.0f, 0.0f};
    geomlib_vec3 v1 = {2.0f, 0.0f, 0.0f};
    geomlib_vec3 v2 = {0.0f, 2.0f, 0.0f};

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 123, 0);

    geomlib_vec3 bary;
    geomlib_vec3 sample = geomlib_sample_triangle_barycentric(v0, v1, v2, &rng, &bary);

    /* Check barycentric coordinates sum to 1 */
    float bary_sum = bary.x + bary.y + bary.z;
    assert(fabsf(bary_sum - 1.0f) < TEST_EPSILON);

    /* Check all barycentric coordinates are non-negative */
    assert(bary.x >= -TEST_EPSILON);
    assert(bary.y >= -TEST_EPSILON);
    assert(bary.z >= -TEST_EPSILON);

    /* Verify point reconstruction using barycentric coordinates */
    geomlib_vec3 reconstructed;
    reconstructed.x = bary.x * v0.x + bary.y * v1.x + bary.z * v2.x;
    reconstructed.y = bary.x * v0.y + bary.y * v1.y + bary.z * v2.y;
    reconstructed.z = bary.x * v0.z + bary.y * v1.z + bary.z * v2.z;

    assert(fabsf(sample.x - reconstructed.x) < TEST_EPSILON);
    assert(fabsf(sample.y - reconstructed.y) < TEST_EPSILON);
    assert(fabsf(sample.z - reconstructed.z) < TEST_EPSILON);

    printf("  ✓ Barycentric coordinates valid (%.4f, %.4f, %.4f sum=%.4f)\n",
           bary.x, bary.y, bary.z, bary_sum);
}

void test_triangle_sample_pdf() {
    printf("Testing triangle sampling PDF...\n");

    /* Unit right triangle has area = 0.5 */
    geomlib_vec3 v0 = {0.0f, 0.0f, 0.0f};
    geomlib_vec3 v1 = {1.0f, 0.0f, 0.0f};
    geomlib_vec3 v2 = {0.0f, 1.0f, 0.0f};

    float pdf = geomlib_triangle_sample_pdf(v0, v1, v2);
    float expected_pdf = 1.0f / 0.5f;  /* 1 / area = 2.0 */

    assert(fabsf(pdf - expected_pdf) < TEST_EPSILON);

    printf("  ✓ PDF correct: %.4f (expected %.4f)\n", pdf, expected_pdf);
}

void test_mesh_sampling_distribution() {
    printf("Testing mesh sampling distribution...\n");

    /* Create a simple mesh with two triangles */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 2);

    geomlib_vec3 vertices[] = {
        {0.0f, 0.0f, 0.0f},
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {1.0f, 1.0f, 0.0f}
    };

    for (int i = 0; i < 4; i++) {
        geomlib_trimesh_add_vertex(&mesh, vertices[i], NULL, NULL);
    }

    /* Triangle 0: vertices 0,1,2 (area = 0.5) */
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    /* Triangle 1: vertices 1,3,2 (area = 0.5) */
    geomlib_trimesh_add_triangle(&mesh, 1, 3, 2);

    geomlib_mesh_sample_distribution dist;
    int result = geomlib_mesh_build_sample_distribution(&mesh, &dist);

    assert(result == GEOMLIB_MESH_SUCCESS);
    assert(dist.triangle_count == 2);
    assert(fabsf(dist.total_area - 1.0f) < TEST_EPSILON);

    /* Each triangle has area 0.5, so equal probability */
    assert(fabsf(dist.triangle_areas[0] - 0.5f) < TEST_EPSILON);
    assert(fabsf(dist.triangle_areas[1] - 0.5f) < TEST_EPSILON);

    /* CDF should be [0.5, 1.0] */
    assert(fabsf(dist.cdf[0] - 0.5f) < TEST_EPSILON);
    assert(fabsf(dist.cdf[1] - 1.0f) < TEST_EPSILON);

    geomlib_mesh_free_sample_distribution(&dist);
    geomlib_trimesh_free(&mesh);

    printf("  ✓ Distribution built correctly (total_area=%.4f)\n", dist.total_area);
}

void test_mesh_surface_sampling() {
    printf("Testing mesh surface sampling...\n");

    /* Create a unit square mesh */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 2);

    geomlib_vec3 vertices[] = {
        {0.0f, 0.0f, 0.0f},
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {1.0f, 1.0f, 0.0f}
    };

    for (int i = 0; i < 4; i++) {
        geomlib_trimesh_add_vertex(&mesh, vertices[i], NULL, NULL);
    }

    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 1, 3, 2);

    geomlib_mesh_sample_distribution dist;
    geomlib_mesh_build_sample_distribution(&mesh, &dist);

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 10000;
    geomlib_vec3 centroid = {0, 0, 0};

    for (int i = 0; i < num_samples; i++) {
        size_t tri_idx;
        geomlib_vec3 sample = geomlib_sample_mesh_surface(&mesh, &dist, &rng, &tri_idx);

        /* Check sample is on the plane z=0 */
        assert(fabsf(sample.z) < TEST_EPSILON);

        /* Check sample is inside unit square */
        assert(sample.x >= -TEST_EPSILON && sample.x <= 1.0f + TEST_EPSILON);
        assert(sample.y >= -TEST_EPSILON && sample.y <= 1.0f + TEST_EPSILON);

        /* Check triangle index is valid */
        assert(tri_idx < mesh.triangle_count);

        centroid.x += sample.x;
        centroid.y += sample.y;
    }

    centroid.x /= num_samples;
    centroid.y /= num_samples;

    /* Expected centroid of unit square is (0.5, 0.5) */
    assert(fabsf(centroid.x - 0.5f) < 0.02f);
    assert(fabsf(centroid.y - 0.5f) < 0.02f);

    geomlib_mesh_free_sample_distribution(&dist);
    geomlib_trimesh_free(&mesh);

    printf("  ✓ Mesh sampling produces uniform distribution (centroid: %.4f, %.4f)\n",
           centroid.x, centroid.y);
}

void test_mc_integration() {
    printf("Testing Monte Carlo integration...\n");

    /* Estimate area of unit circle using Monte Carlo integration */
    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 100000;
    float *values = (float *)malloc(sizeof(float) * num_samples);
    geomlib_vec3 *samples = (geomlib_vec3 *)malloc(sizeof(geomlib_vec3) * num_samples);

    /* Sample points in unit square [-1,1]^2 */
    for (int i = 0; i < num_samples; i++) {
        samples[i].x = geomlib_pcg32_range_float(&rng, -1.0f, 1.0f);
        samples[i].y = geomlib_pcg32_range_float(&rng, -1.0f, 1.0f);
        samples[i].z = 0.0f;

        /* Indicator function: 1 if inside unit circle, 0 otherwise */
        float dist_sq = samples[i].x * samples[i].x + samples[i].y * samples[i].y;
        values[i] = (dist_sq <= 1.0f) ? 1.0f : 0.0f;
    }

    /* Area of sampling domain (square) is 4.0 */
    double estimate = geomlib_mc_estimate_surface_integral(
        samples, values, NULL, num_samples, 4.0f
    );

    /* Expected value is π ≈ 3.14159 */
    double expected = M_PI;
    double error = fabs(estimate - expected);

    /* With 100k samples, error should be small */
    assert(error < 0.05);

    free(values);
    free(samples);

    printf("  ✓ Monte Carlo integration estimate: %.5f (expected %.5f, error %.5f)\n",
           estimate, expected, error);
}

void test_closest_point_triangle() {
    printf("Testing closest point on triangle...\n");

    geomlib_vec3 v0 = {0.0f, 0.0f, 0.0f};
    geomlib_vec3 v1 = {1.0f, 0.0f, 0.0f};
    geomlib_vec3 v2 = {0.0f, 1.0f, 0.0f};

    /* Test 1: Point inside triangle */
    geomlib_vec3 query1 = {0.25f, 0.25f, 0.0f};
    geomlib_vec3 closest1 = geomlib_closest_point_triangle(query1, v0, v1, v2);
    assert(fabsf(closest1.x - 0.25f) < TEST_EPSILON);
    assert(fabsf(closest1.y - 0.25f) < TEST_EPSILON);
    assert(fabsf(closest1.z) < TEST_EPSILON);

    /* Test 2: Point above triangle */
    geomlib_vec3 query2 = {0.25f, 0.25f, 1.0f};
    geomlib_vec3 closest2 = geomlib_closest_point_triangle(query2, v0, v1, v2);
    assert(fabsf(closest2.x - 0.25f) < TEST_EPSILON);
    assert(fabsf(closest2.y - 0.25f) < TEST_EPSILON);
    assert(fabsf(closest2.z) < TEST_EPSILON);

    /* Test 3: Point outside, closest to vertex */
    geomlib_vec3 query3 = {-1.0f, -1.0f, 0.0f};
    geomlib_vec3 closest3 = geomlib_closest_point_triangle(query3, v0, v1, v2);
    assert(fabsf(closest3.x) < TEST_EPSILON);
    assert(fabsf(closest3.y) < TEST_EPSILON);
    assert(fabsf(closest3.z) < TEST_EPSILON);

    printf("  ✓ Closest point queries work correctly\n");
}

void test_sphere_sampling() {
    printf("Testing sphere surface sampling...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 10000;
    geomlib_vec3 centroid = {0, 0, 0};

    for (int i = 0; i < num_samples; i++) {
        geomlib_vec3 sample = geomlib_sample_sphere_uniform(&rng);

        /* Check sample is on unit sphere */
        float length = geomlib_vec3_length(sample);
        assert(fabsf(length - 1.0f) < TEST_EPSILON);

        centroid.x += sample.x;
        centroid.y += sample.y;
        centroid.z += sample.z;
    }

    centroid.x /= num_samples;
    centroid.y /= num_samples;
    centroid.z /= num_samples;

    /* Centroid should be near origin for uniform sampling */
    assert(fabsf(centroid.x) < 0.02f);
    assert(fabsf(centroid.y) < 0.02f);
    assert(fabsf(centroid.z) < 0.02f);

    printf("  ✓ Uniform sphere sampling correct (centroid near origin: %.4f, %.4f, %.4f)\n",
           centroid.x, centroid.y, centroid.z);
}

void test_sphere_volume_sampling() {
    printf("Testing sphere volume sampling...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 10000;
    int inside_count = 0;

    for (int i = 0; i < num_samples; i++) {
        geomlib_vec3 sample = geomlib_sample_sphere_volume(&rng);

        /* Check sample is inside unit sphere */
        float length_sq = sample.x * sample.x + sample.y * sample.y + sample.z * sample.z;
        assert(length_sq <= 1.0f + TEST_EPSILON);

        if (length_sq <= 1.0f) {
            inside_count++;
        }
    }

    /* All samples should be inside */
    assert(inside_count == num_samples);

    printf("  ✓ All %d samples inside unit sphere\n", inside_count);
}

void test_hemisphere_sampling() {
    printf("Testing hemisphere sampling...\n");

    geomlib_vec3 normal = {0.0f, 0.0f, 1.0f};
    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 1000;
    int above_plane_count = 0;

    for (int i = 0; i < num_samples; i++) {
        geomlib_vec3 dir = geomlib_sample_hemisphere_uniform(normal, &rng);

        /* Check direction is normalized */
        float length = geomlib_vec3_length(dir);
        assert(fabsf(length - 1.0f) < TEST_EPSILON);

        /* Check direction is in upper hemisphere (dot with normal >= 0) */
        float dot = geomlib_vec3_dot(dir, normal);
        assert(dot >= -TEST_EPSILON);

        if (dot > 0.0f) {
            above_plane_count++;
        }
    }

    /* Most samples should be strictly above the plane */
    assert(above_plane_count > num_samples * 0.95);

    printf("  ✓ Hemisphere sampling correct (%d/%d above plane)\n", above_plane_count, num_samples);
}

void test_cosine_hemisphere_sampling() {
    printf("Testing cosine-weighted hemisphere sampling...\n");

    geomlib_vec3 normal = {0.0f, 1.0f, 0.0f};
    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int num_samples = 10000;
    double sum_cos = 0.0;

    for (int i = 0; i < num_samples; i++) {
        geomlib_vec3 dir = geomlib_sample_hemisphere_cosine(normal, &rng);

        /* Check direction is normalized */
        float length = geomlib_vec3_length(dir);
        assert(fabsf(length - 1.0f) < TEST_EPSILON);

        /* Check direction is in hemisphere */
        float cos_theta = geomlib_vec3_dot(dir, normal);
        assert(cos_theta >= -TEST_EPSILON);

        sum_cos += cos_theta;
    }

    double avg_cos = sum_cos / num_samples;

    /* For cosine-weighted sampling, average cos(theta) should be around 2/3 */
    /* This comes from: ∫∫ cos(θ) * cos(θ) sin(θ) dθ dφ / ∫∫ cos(θ) sin(θ) dθ dφ = 2/3 */
    assert(avg_cos > 0.6 && avg_cos < 0.7);

    printf("  ✓ Cosine hemisphere sampling correct (avg cos = %.4f, expected ≈0.67)\n", avg_cos);
}

int main() {
    printf("================================================\n");
    printf("geomlib Monte Carlo Geometry Processing Tests\n");
    printf("================================================\n\n");

    test_triangle_sampling_uniform();
    test_triangle_sampling_barycentric();
    test_triangle_sample_pdf();
    test_mesh_sampling_distribution();
    test_mesh_surface_sampling();
    test_mc_integration();
    test_closest_point_triangle();
    test_sphere_sampling();
    test_sphere_volume_sampling();
    test_hemisphere_sampling();
    test_cosine_hemisphere_sampling();

    printf("\n================================================\n");
    printf("All Monte Carlo tests passed!\n");
    printf("================================================\n");

    return 0;
}
