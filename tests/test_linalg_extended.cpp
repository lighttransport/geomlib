/*
 * geomlib - Extended Linear Algebra Tests (Simplified)
 */

#include <cstdio>
#include <cmath>
#include "geomlib/geomlib_linalg.h"

static int test_count = 0;
static int pass_count = 0;

#define EPSILON 1e-5f

static bool float_eq(float a, float b, float eps = EPSILON) {
    return std::fabs(a - b) < eps;
}

static void test(bool condition, const char* name) {
    test_count++;
    if (condition) {
        pass_count++;
        std::printf("[PASS] %s\n", name);
    } else {
        std::printf("[FAIL] %s\n", name);
    }
}

static void test_vec2d_operations() {
    geomlib_vec2d a = geomlib_vec2d_create(3.0, 4.0);
    geomlib_vec2d b = geomlib_vec2d_create(1.0, 2.0);

    geomlib_vec2d sum = geomlib_vec2d_add(a, b);
    test(float_eq(sum.x, 4.0) && float_eq(sum.y, 6.0), "vec2d_add");

    double length = geomlib_vec2d_length(a);
    test(float_eq(length, 5.0, 1e-6), "vec2d_length");

    geomlib_vec2d norm = geomlib_vec2d_normalize(a);
    test(float_eq(norm.x, 0.6, 1e-6) && float_eq(norm.y, 0.8, 1e-6), "vec2d_normalize");
}

static void test_vec3_edge_cases() {
    geomlib_vec3 zero = geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    test(float_eq(geomlib_vec3_length(zero), 0.0f), "vec3_zero_length");

    geomlib_vec3 v1 = geomlib_vec3_create(1.0f, 0.0f, 0.0f);
    geomlib_vec3 v2 = geomlib_vec3_create(2.0f, 0.0f, 0.0f);
    geomlib_vec3 cross = geomlib_vec3_cross(v1, v2);
    test(float_eq(cross.x, 0.0f) && float_eq(cross.y, 0.0f) && float_eq(cross.z, 0.0f),
         "vec3_parallel_cross_zero");
}

static void test_vec4_operations() {
    geomlib_vec4 a = geomlib_vec4_create(1.0f, 2.0f, 3.0f, 1.0f);
    geomlib_vec4 b = geomlib_vec4_create(4.0f, 5.0f, 6.0f, 1.0f);

    geomlib_vec4 sum = geomlib_vec4_add(a, b);
    test(float_eq(sum.x, 5.0f) && float_eq(sum.y, 7.0f), "vec4_add");

    geomlib_vec4 scaled = geomlib_vec4_mul(a, 2.0f);
    test(float_eq(scaled.x, 2.0f) && float_eq(scaled.y, 4.0f), "vec4_mul");

    float dot = geomlib_vec4_dot(a, b);
    test(float_eq(dot, 33.0f), "vec4_dot");
}

static void test_mat4_transformations() {
    geomlib_mat4 translate = geomlib_mat4_translate(1.0f, 2.0f, 3.0f);
    geomlib_mat4 scale = geomlib_mat4_scale_xyz(2.0f, 2.0f, 2.0f);

    geomlib_vec4 point = geomlib_vec4_create(1.0f, 1.0f, 1.0f, 1.0f);
    geomlib_vec4 scaled = geomlib_mat4_mul_vec4(scale, point);
    test(float_eq(scaled.x, 2.0f) && float_eq(scaled.y, 2.0f), "mat4_scale_xyz");
}

static void test_mat4_inverse() {
    geomlib_mat4 m = geomlib_mat4_translate(5.0f, 3.0f, 2.0f);
    geomlib_mat4 inv;
    int result = geomlib_mat4_inverse(m, &inv);
    test(result == 1, "mat4_inverse_succeeds");
}

static void test_vec3_lerp() {
    geomlib_vec3 a = geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    geomlib_vec3 b = geomlib_vec3_create(10.0f, 10.0f, 10.0f);

    geomlib_vec3 result = geomlib_vec3_lerp(a, b, 0.0f);
    test(float_eq(result.x, 0.0f), "vec3_lerp_t0");

    result = geomlib_vec3_lerp(a, b, 1.0f);
    test(float_eq(result.x, 10.0f), "vec3_lerp_t1");

    result = geomlib_vec3_lerp(a, b, 0.5f);
    test(float_eq(result.x, 5.0f), "vec3_lerp_t05");
}

int main() {
    std::printf("=== Extended Linear Algebra Tests ===\n\n");

    test_vec2d_operations();
    test_vec3_edge_cases();
    test_vec4_operations();
    test_mat4_transformations();
    test_mat4_inverse();
    test_vec3_lerp();

    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);

    return (pass_count == test_count) ? 0 : 1;
}
