/*
 * geomlib - Linear Algebra Tests
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

static void test_vec2() {
    geomlib_vec2 a = geomlib_vec2_create(1.0f, 2.0f);
    geomlib_vec2 b = geomlib_vec2_create(3.0f, 4.0f);
    
    test(float_eq(a.x, 1.0f) && float_eq(a.y, 2.0f), "vec2_create");
    
    geomlib_vec2 sum = geomlib_vec2_add(a, b);
    test(float_eq(sum.x, 4.0f) && float_eq(sum.y, 6.0f), "vec2_add");
    
    geomlib_vec2 diff = geomlib_vec2_sub(a, b);
    test(float_eq(diff.x, -2.0f) && float_eq(diff.y, -2.0f), "vec2_sub");
    
    geomlib_vec2 scaled = geomlib_vec2_mul(a, 2.0f);
    test(float_eq(scaled.x, 2.0f) && float_eq(scaled.y, 4.0f), "vec2_mul");
    
    float dot = geomlib_vec2_dot(a, b);
    test(float_eq(dot, 11.0f), "vec2_dot");  // 1*3 + 2*4 = 11
    
    float len = geomlib_vec2_length(geomlib_vec2_create(3.0f, 4.0f));
    test(float_eq(len, 5.0f), "vec2_length");  // 3-4-5 triangle
    
    geomlib_vec2 norm = geomlib_vec2_normalize(geomlib_vec2_create(3.0f, 4.0f));
    test(float_eq(norm.x, 0.6f) && float_eq(norm.y, 0.8f), "vec2_normalize");
}

static void test_vec3() {
    geomlib_vec3 a = geomlib_vec3_create(1.0f, 2.0f, 3.0f);
    geomlib_vec3 b = geomlib_vec3_create(4.0f, 5.0f, 6.0f);
    
    test(float_eq(a.x, 1.0f) && float_eq(a.y, 2.0f) && float_eq(a.z, 3.0f), "vec3_create");
    
    geomlib_vec3 sum = geomlib_vec3_add(a, b);
    test(float_eq(sum.x, 5.0f) && float_eq(sum.y, 7.0f) && float_eq(sum.z, 9.0f), "vec3_add");
    
    float dot = geomlib_vec3_dot(a, b);
    test(float_eq(dot, 32.0f), "vec3_dot");  // 1*4 + 2*5 + 3*6 = 32
    
    geomlib_vec3 i = geomlib_vec3_create(1.0f, 0.0f, 0.0f);
    geomlib_vec3 j = geomlib_vec3_create(0.0f, 1.0f, 0.0f);
    geomlib_vec3 k = geomlib_vec3_cross(i, j);
    test(float_eq(k.x, 0.0f) && float_eq(k.y, 0.0f) && float_eq(k.z, 1.0f), "vec3_cross");
    
    geomlib_vec3 lerp = geomlib_vec3_lerp(a, b, 0.5f);
    test(float_eq(lerp.x, 2.5f) && float_eq(lerp.y, 3.5f) && float_eq(lerp.z, 4.5f), "vec3_lerp");
}

static void test_vec4() {
    geomlib_vec4 a = geomlib_vec4_create(1.0f, 2.0f, 3.0f, 4.0f);
    geomlib_vec4 b = geomlib_vec4_create(5.0f, 6.0f, 7.0f, 8.0f);
    
    float dot = geomlib_vec4_dot(a, b);
    test(float_eq(dot, 70.0f), "vec4_dot");  // 1*5 + 2*6 + 3*7 + 4*8 = 70
    
    float len = geomlib_vec4_length(geomlib_vec4_create(1.0f, 0.0f, 0.0f, 0.0f));
    test(float_eq(len, 1.0f), "vec4_length");
}

static void test_mat3() {
    geomlib_mat3 I = geomlib_mat3_identity();
    test(float_eq(I.m[0], 1.0f) && float_eq(I.m[4], 1.0f) && float_eq(I.m[8], 1.0f), "mat3_identity");
    
    geomlib_mat3 A = geomlib_mat3_identity();
    geomlib_mat3 B = geomlib_mat3_mul(A, I);
    test(float_eq(B.m[0], 1.0f) && float_eq(B.m[4], 1.0f) && float_eq(B.m[8], 1.0f), "mat3_mul_identity");
    
    geomlib_vec3 v = geomlib_vec3_create(1.0f, 2.0f, 3.0f);
    geomlib_vec3 result = geomlib_mat3_mul_vec3(I, v);
    test(float_eq(result.x, 1.0f) && float_eq(result.y, 2.0f) && float_eq(result.z, 3.0f), "mat3_mul_vec3_identity");
    
    float det = geomlib_mat3_determinant(I);
    test(float_eq(det, 1.0f), "mat3_determinant_identity");
    
    geomlib_mat3 inv;
    int success = geomlib_mat3_inverse(I, &inv);
    test(success && float_eq(inv.m[0], 1.0f), "mat3_inverse_identity");
}

static void test_mat4() {
    geomlib_mat4 I = geomlib_mat4_identity();
    test(float_eq(I.m[0], 1.0f) && float_eq(I.m[5], 1.0f) && float_eq(I.m[10], 1.0f) && float_eq(I.m[15], 1.0f), 
         "mat4_identity");
    
    geomlib_mat4 T = geomlib_mat4_translate(1.0f, 2.0f, 3.0f);
    geomlib_vec4 p = geomlib_vec4_create(0.0f, 0.0f, 0.0f, 1.0f);
    geomlib_vec4 result = geomlib_mat4_mul_vec4(T, p);
    test(float_eq(result.x, 1.0f) && float_eq(result.y, 2.0f) && float_eq(result.z, 3.0f), "mat4_translate");
    
    geomlib_mat4 S = geomlib_mat4_scale_xyz(2.0f, 3.0f, 4.0f);
    p = geomlib_vec4_create(1.0f, 1.0f, 1.0f, 1.0f);
    result = geomlib_mat4_mul_vec4(S, p);
    test(float_eq(result.x, 2.0f) && float_eq(result.y, 3.0f) && float_eq(result.z, 4.0f), "mat4_scale");
    
    float det = geomlib_mat4_determinant(I);
    test(float_eq(det, 1.0f), "mat4_determinant_identity");
    
    geomlib_mat4 inv;
    int success = geomlib_mat4_inverse(I, &inv);
    test(success && float_eq(inv.m[0], 1.0f), "mat4_inverse_identity");
    
    // Test rotation
    geomlib_mat4 Rx = geomlib_mat4_rotate_x(0.0f);
    geomlib_vec4 v = geomlib_vec4_create(0.0f, 1.0f, 0.0f, 0.0f);
    result = geomlib_mat4_mul_vec4(Rx, v);
    test(float_eq(result.y, 1.0f), "mat4_rotate_x_zero");
}

static void test_mat4_inverse() {
    // Test a more complex matrix inverse
    geomlib_mat4 A;
    for (int i = 0; i < 16; i++) A.m[i] = 0.0f;
    A.m[0] = 1.0f; A.m[5] = 2.0f; A.m[10] = 3.0f; A.m[15] = 1.0f;
    A.m[3] = 1.0f;  // Translation component
    
    geomlib_mat4 inv;
    int success = geomlib_mat4_inverse(A, &inv);
    
    // A * inv should be identity
    geomlib_mat4 product = geomlib_mat4_mul(A, inv);
    
    bool is_identity = true;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            float expected = (i == j) ? 1.0f : 0.0f;
            if (!float_eq(product.m[j * 4 + i], expected, 1e-4f)) {
                is_identity = false;
            }
        }
    }
    
    test(success && is_identity, "mat4_inverse_product_is_identity");
}

int main() {
    std::printf("=== Linear Algebra Tests ===\n\n");
    
    test_vec2();
    test_vec3();
    test_vec4();
    test_mat3();
    test_mat4();
    test_mat4_inverse();
    
    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);
    
    return (pass_count == test_count) ? 0 : 1;
}
