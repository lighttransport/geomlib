/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Linear Algebra Implementation
 */

#include "geomlib/geomlib_allocator.h"
#include "geomlib/geomlib_c11.h"
#include "geomlib/geomlib_linalg.h"

/* Column-major matrix access helper - element at row r, col c */
/* mat3: index = c*3 + r */
/* mat4: index = c*4 + r */

/* ============================================================================
 * Vector2 Operations (float)
 * ============================================================================ */

geomlib_vec2 geomlib_vec2_create(float x, float y) {
    geomlib_vec2 v;
    v.x = x;
    v.y = y;
    return v;
}

geomlib_vec2 geomlib_vec2_add(geomlib_vec2 va, geomlib_vec2 vb) {
    geomlib_vec2 v;
    v.x = va.x + vb.x;
    v.y = va.y + vb.y;
    return v;
}

geomlib_vec2 geomlib_vec2_sub(geomlib_vec2 va, geomlib_vec2 vb) {
    geomlib_vec2 v;
    v.x = va.x - vb.x;
    v.y = va.y - vb.y;
    return v;
}

geomlib_vec2 geomlib_vec2_mul(geomlib_vec2 v, float s) {
    geomlib_vec2 result;
    result.x = v.x * s;
    result.y = v.y * s;
    return result;
}

geomlib_vec2 geomlib_vec2_div(geomlib_vec2 v, float s) {
    geomlib_vec2 result;
    float inv_s = 1.0f / s;
    result.x = v.x * inv_s;
    result.y = v.y * inv_s;
    return result;
}

float geomlib_vec2_dot(geomlib_vec2 va, geomlib_vec2 vb) {
    return va.x * vb.x + va.y * vb.y;
}

float geomlib_vec2_length_squared(geomlib_vec2 v) {
    return v.x * v.x + v.y * v.y;
}

float geomlib_vec2_length(geomlib_vec2 v) {
    return sqrtf(geomlib_vec2_length_squared(v));
}

geomlib_vec2 geomlib_vec2_normalize(geomlib_vec2 v) {
    float len = geomlib_vec2_length(v);
    if (len > 0.0f) {
        return geomlib_vec2_div(v, len);
    }
    return v;
}

geomlib_vec2 geomlib_vec2_negate(geomlib_vec2 v) {
    geomlib_vec2 result;
    result.x = -v.x;
    result.y = -v.y;
    return result;
}

/* ============================================================================
 * Vector3 Operations (float)
 * ============================================================================ */

geomlib_vec3 geomlib_vec3_create(float x, float y, float z) {
    geomlib_vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

geomlib_vec3 geomlib_vec3_add(geomlib_vec3 va, geomlib_vec3 vb) {
    geomlib_vec3 v;
    v.x = va.x + vb.x;
    v.y = va.y + vb.y;
    v.z = va.z + vb.z;
    return v;
}

geomlib_vec3 geomlib_vec3_sub(geomlib_vec3 va, geomlib_vec3 vb) {
    geomlib_vec3 v;
    v.x = va.x - vb.x;
    v.y = va.y - vb.y;
    v.z = va.z - vb.z;
    return v;
}

geomlib_vec3 geomlib_vec3_mul(geomlib_vec3 v, float s) {
    geomlib_vec3 result;
    result.x = v.x * s;
    result.y = v.y * s;
    result.z = v.z * s;
    return result;
}

geomlib_vec3 geomlib_vec3_div(geomlib_vec3 v, float s) {
    geomlib_vec3 result;
    float inv_s = 1.0f / s;
    result.x = v.x * inv_s;
    result.y = v.y * inv_s;
    result.z = v.z * inv_s;
    return result;
}

float geomlib_vec3_dot(geomlib_vec3 va, geomlib_vec3 vb) {
    return va.x * vb.x + va.y * vb.y + va.z * vb.z;
}

geomlib_vec3 geomlib_vec3_cross(geomlib_vec3 va, geomlib_vec3 vb) {
    geomlib_vec3 v;
    v.x = va.y * vb.z - va.z * vb.y;
    v.y = va.z * vb.x - va.x * vb.z;
    v.z = va.x * vb.y - va.y * vb.x;
    return v;
}

float geomlib_vec3_length_squared(geomlib_vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

float geomlib_vec3_length(geomlib_vec3 v) {
    return sqrtf(geomlib_vec3_length_squared(v));
}

geomlib_vec3 geomlib_vec3_normalize(geomlib_vec3 v) {
    float len = geomlib_vec3_length(v);
    if (len > 0.0f) {
        return geomlib_vec3_div(v, len);
    }
    return v;
}

geomlib_vec3 geomlib_vec3_negate(geomlib_vec3 v) {
    geomlib_vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;
    return result;
}

geomlib_vec3 geomlib_vec3_lerp(geomlib_vec3 va, geomlib_vec3 vb, float t) {
    geomlib_vec3 v;
    v.x = va.x + (vb.x - va.x) * t;
    v.y = va.y + (vb.y - va.y) * t;
    v.z = va.z + (vb.z - va.z) * t;
    return v;
}

/* ============================================================================
 * Vector4 Operations (float)
 * ============================================================================ */

geomlib_vec4 geomlib_vec4_create(float x, float y, float z, float w) {
    geomlib_vec4 v;
    v.x = x;
    v.y = y;
    v.z = z;
    v.w = w;
    return v;
}

geomlib_vec4 geomlib_vec4_add(geomlib_vec4 va, geomlib_vec4 vb) {
    geomlib_vec4 v;
    v.x = va.x + vb.x;
    v.y = va.y + vb.y;
    v.z = va.z + vb.z;
    v.w = va.w + vb.w;
    return v;
}

geomlib_vec4 geomlib_vec4_sub(geomlib_vec4 va, geomlib_vec4 vb) {
    geomlib_vec4 v;
    v.x = va.x - vb.x;
    v.y = va.y - vb.y;
    v.z = va.z - vb.z;
    v.w = va.w - vb.w;
    return v;
}

geomlib_vec4 geomlib_vec4_mul(geomlib_vec4 v, float s) {
    geomlib_vec4 result;
    result.x = v.x * s;
    result.y = v.y * s;
    result.z = v.z * s;
    result.w = v.w * s;
    return result;
}

geomlib_vec4 geomlib_vec4_div(geomlib_vec4 v, float s) {
    geomlib_vec4 result;
    float inv_s = 1.0f / s;
    result.x = v.x * inv_s;
    result.y = v.y * inv_s;
    result.z = v.z * inv_s;
    result.w = v.w * inv_s;
    return result;
}

float geomlib_vec4_dot(geomlib_vec4 va, geomlib_vec4 vb) {
    return va.x * vb.x + va.y * vb.y + va.z * vb.z + va.w * vb.w;
}

float geomlib_vec4_length_squared(geomlib_vec4 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
}

float geomlib_vec4_length(geomlib_vec4 v) {
    return sqrtf(geomlib_vec4_length_squared(v));
}

geomlib_vec4 geomlib_vec4_normalize(geomlib_vec4 v) {
    float len = geomlib_vec4_length(v);
    if (len > 0.0f) {
        return geomlib_vec4_div(v, len);
    }
    return v;
}

geomlib_vec4 geomlib_vec4_negate(geomlib_vec4 v) {
    geomlib_vec4 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;
    result.w = -v.w;
    return result;
}

/* ============================================================================
 * Matrix3 Operations (float) - Column-major
 * ============================================================================ */

geomlib_mat3 geomlib_mat3_identity(void) {
    geomlib_mat3 mat;
    geomlib_memset(&mat, 0, sizeof(mat));
    mat.m[0] = 1.0f;  /* (0,0) */
    mat.m[4] = 1.0f;  /* (1,1) */
    mat.m[8] = 1.0f;  /* (2,2) */
    return mat;
}

geomlib_mat3 geomlib_mat3_zero(void) {
    geomlib_mat3 mat;
    geomlib_memset(&mat, 0, sizeof(mat));
    return mat;
}

geomlib_mat3 geomlib_mat3_add(geomlib_mat3 ma, geomlib_mat3 mb) {
    geomlib_mat3 mat;
    for (int idx = 0; idx < 9; idx++) {
        mat.m[idx] = ma.m[idx] + mb.m[idx];
    }
    return mat;
}

geomlib_mat3 geomlib_mat3_sub(geomlib_mat3 ma, geomlib_mat3 mb) {
    geomlib_mat3 mat;
    for (int idx = 0; idx < 9; idx++) {
        mat.m[idx] = ma.m[idx] - mb.m[idx];
    }
    return mat;
}

geomlib_mat3 geomlib_mat3_mul(geomlib_mat3 ma, geomlib_mat3 mb) {
    geomlib_mat3 mat;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < 3; row++) {
            mat.m[col * 3 + row] = 
                ma.m[0 * 3 + row] * mb.m[col * 3 + 0] +
                ma.m[1 * 3 + row] * mb.m[col * 3 + 1] +
                ma.m[2 * 3 + row] * mb.m[col * 3 + 2];
        }
    }
    return mat;
}

geomlib_mat3 geomlib_mat3_scale(geomlib_mat3 mat, float s) {
    geomlib_mat3 result;
    for (int idx = 0; idx < 9; idx++) {
        result.m[idx] = mat.m[idx] * s;
    }
    return result;
}

geomlib_vec3 geomlib_mat3_mul_vec3(geomlib_mat3 mat, geomlib_vec3 v) {
    geomlib_vec3 result;
    result.x = mat.m[0 * 3 + 0] * v.x + mat.m[1 * 3 + 0] * v.y + mat.m[2 * 3 + 0] * v.z;
    result.y = mat.m[0 * 3 + 1] * v.x + mat.m[1 * 3 + 1] * v.y + mat.m[2 * 3 + 1] * v.z;
    result.z = mat.m[0 * 3 + 2] * v.x + mat.m[1 * 3 + 2] * v.y + mat.m[2 * 3 + 2] * v.z;
    return result;
}

geomlib_mat3 geomlib_mat3_transpose(geomlib_mat3 mat) {
    geomlib_mat3 result;
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            result.m[col * 3 + row] = mat.m[row * 3 + col];
        }
    }
    return result;
}

float geomlib_mat3_determinant(geomlib_mat3 mat) {
    float m00 = mat.m[0 * 3 + 0], m01 = mat.m[1 * 3 + 0], m02 = mat.m[2 * 3 + 0];
    float m10 = mat.m[0 * 3 + 1], m11 = mat.m[1 * 3 + 1], m12 = mat.m[2 * 3 + 1];
    float m20 = mat.m[0 * 3 + 2], m21 = mat.m[1 * 3 + 2], m22 = mat.m[2 * 3 + 2];
    
    return m00 * (m11 * m22 - m21 * m12) -
           m01 * (m10 * m22 - m12 * m20) +
           m02 * (m10 * m21 - m11 * m20);
}

int geomlib_mat3_inverse(geomlib_mat3 mat, geomlib_mat3 *out) {
    float det = geomlib_mat3_determinant(mat);
    if (fabsf(det) < 1e-10f) {
        return 0; /* Singular matrix */
    }
    
    float inv_det = 1.0f / det;
    
    float m00 = mat.m[0 * 3 + 0], m01 = mat.m[1 * 3 + 0], m02 = mat.m[2 * 3 + 0];
    float m10 = mat.m[0 * 3 + 1], m11 = mat.m[1 * 3 + 1], m12 = mat.m[2 * 3 + 1];
    float m20 = mat.m[0 * 3 + 2], m21 = mat.m[1 * 3 + 2], m22 = mat.m[2 * 3 + 2];
    
    out->m[0 * 3 + 0] = (m11 * m22 - m21 * m12) * inv_det;
    out->m[1 * 3 + 0] = (m02 * m21 - m01 * m22) * inv_det;
    out->m[2 * 3 + 0] = (m01 * m12 - m02 * m11) * inv_det;
    out->m[0 * 3 + 1] = (m12 * m20 - m10 * m22) * inv_det;
    out->m[1 * 3 + 1] = (m00 * m22 - m02 * m20) * inv_det;
    out->m[2 * 3 + 1] = (m10 * m02 - m00 * m12) * inv_det;
    out->m[0 * 3 + 2] = (m10 * m21 - m20 * m11) * inv_det;
    out->m[1 * 3 + 2] = (m20 * m01 - m00 * m21) * inv_det;
    out->m[2 * 3 + 2] = (m00 * m11 - m10 * m01) * inv_det;
    
    return 1;
}

/* ============================================================================
 * Matrix4 Operations (float) - Column-major
 * ============================================================================ */

geomlib_mat4 geomlib_mat4_identity(void) {
    geomlib_mat4 mat;
    geomlib_memset(&mat, 0, sizeof(mat));
    mat.m[0] = 1.0f;   /* (0,0) */
    mat.m[5] = 1.0f;   /* (1,1) */
    mat.m[10] = 1.0f;  /* (2,2) */
    mat.m[15] = 1.0f;  /* (3,3) */
    return mat;
}

geomlib_mat4 geomlib_mat4_zero(void) {
    geomlib_mat4 mat;
    geomlib_memset(&mat, 0, sizeof(mat));
    return mat;
}

geomlib_mat4 geomlib_mat4_add(geomlib_mat4 ma, geomlib_mat4 mb) {
    geomlib_mat4 mat;
    for (int idx = 0; idx < 16; idx++) {
        mat.m[idx] = ma.m[idx] + mb.m[idx];
    }
    return mat;
}

geomlib_mat4 geomlib_mat4_sub(geomlib_mat4 ma, geomlib_mat4 mb) {
    geomlib_mat4 mat;
    for (int idx = 0; idx < 16; idx++) {
        mat.m[idx] = ma.m[idx] - mb.m[idx];
    }
    return mat;
}

geomlib_mat4 geomlib_mat4_mul(geomlib_mat4 ma, geomlib_mat4 mb) {
    geomlib_mat4 mat;
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            mat.m[col * 4 + row] = 
                ma.m[0 * 4 + row] * mb.m[col * 4 + 0] +
                ma.m[1 * 4 + row] * mb.m[col * 4 + 1] +
                ma.m[2 * 4 + row] * mb.m[col * 4 + 2] +
                ma.m[3 * 4 + row] * mb.m[col * 4 + 3];
        }
    }
    return mat;
}

geomlib_mat4 geomlib_mat4_scale(geomlib_mat4 mat, float s) {
    geomlib_mat4 result;
    for (int idx = 0; idx < 16; idx++) {
        result.m[idx] = mat.m[idx] * s;
    }
    return result;
}

geomlib_vec4 geomlib_mat4_mul_vec4(geomlib_mat4 mat, geomlib_vec4 v) {
    geomlib_vec4 result;
    result.x = mat.m[0 * 4 + 0] * v.x + mat.m[1 * 4 + 0] * v.y + mat.m[2 * 4 + 0] * v.z + mat.m[3 * 4 + 0] * v.w;
    result.y = mat.m[0 * 4 + 1] * v.x + mat.m[1 * 4 + 1] * v.y + mat.m[2 * 4 + 1] * v.z + mat.m[3 * 4 + 1] * v.w;
    result.z = mat.m[0 * 4 + 2] * v.x + mat.m[1 * 4 + 2] * v.y + mat.m[2 * 4 + 2] * v.z + mat.m[3 * 4 + 2] * v.w;
    result.w = mat.m[0 * 4 + 3] * v.x + mat.m[1 * 4 + 3] * v.y + mat.m[2 * 4 + 3] * v.z + mat.m[3 * 4 + 3] * v.w;
    return result;
}

geomlib_mat4 geomlib_mat4_transpose(geomlib_mat4 mat) {
    geomlib_mat4 result;
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            result.m[col * 4 + row] = mat.m[row * 4 + col];
        }
    }
    return result;
}

float geomlib_mat4_determinant(geomlib_mat4 mat) {
    float m00 = mat.m[0 * 4 + 0], m01 = mat.m[1 * 4 + 0], m02 = mat.m[2 * 4 + 0], m03 = mat.m[3 * 4 + 0];
    float m10 = mat.m[0 * 4 + 1], m11 = mat.m[1 * 4 + 1], m12 = mat.m[2 * 4 + 1], m13 = mat.m[3 * 4 + 1];
    float m20 = mat.m[0 * 4 + 2], m21 = mat.m[1 * 4 + 2], m22 = mat.m[2 * 4 + 2], m23 = mat.m[3 * 4 + 2];
    float m30 = mat.m[0 * 4 + 3], m31 = mat.m[1 * 4 + 3], m32 = mat.m[2 * 4 + 3], m33 = mat.m[3 * 4 + 3];
    
    float kp_lo = m22 * m33 - m23 * m32;
    float jp_ln = m21 * m33 - m23 * m31;
    float jo_kn = m21 * m32 - m22 * m31;
    float ip_lm = m20 * m33 - m23 * m30;
    float io_km = m20 * m32 - m22 * m30;
    float in_jm = m20 * m31 - m21 * m30;
    
    return m00 * (m11 * kp_lo - m12 * jp_ln + m13 * jo_kn) -
           m01 * (m10 * kp_lo - m12 * ip_lm + m13 * io_km) +
           m02 * (m10 * jp_ln - m11 * ip_lm + m13 * in_jm) -
           m03 * (m10 * jo_kn - m11 * io_km + m12 * in_jm);
}

int geomlib_mat4_inverse(geomlib_mat4 mat, geomlib_mat4 *out) {
    float m00 = mat.m[0 * 4 + 0], m01 = mat.m[1 * 4 + 0], m02 = mat.m[2 * 4 + 0], m03 = mat.m[3 * 4 + 0];
    float m10 = mat.m[0 * 4 + 1], m11 = mat.m[1 * 4 + 1], m12 = mat.m[2 * 4 + 1], m13 = mat.m[3 * 4 + 1];
    float m20 = mat.m[0 * 4 + 2], m21 = mat.m[1 * 4 + 2], m22 = mat.m[2 * 4 + 2], m23 = mat.m[3 * 4 + 2];
    float m30 = mat.m[0 * 4 + 3], m31 = mat.m[1 * 4 + 3], m32 = mat.m[2 * 4 + 3], m33 = mat.m[3 * 4 + 3];
    
    float kp_lo = m22 * m33 - m23 * m32;
    float jp_ln = m21 * m33 - m23 * m31;
    float jo_kn = m21 * m32 - m22 * m31;
    float ip_lm = m20 * m33 - m23 * m30;
    float io_km = m20 * m32 - m22 * m30;
    float in_jm = m20 * m31 - m21 * m30;
    
    float a11 = m11 * kp_lo - m12 * jp_ln + m13 * jo_kn;
    float a12 = -(m10 * kp_lo - m12 * ip_lm + m13 * io_km);
    float a13 = m10 * jp_ln - m11 * ip_lm + m13 * in_jm;
    float a14 = -(m10 * jo_kn - m11 * io_km + m12 * in_jm);
    
    float det = m00 * a11 + m01 * a12 + m02 * a13 + m03 * a14;
    
    if (fabsf(det) < 1e-10f) {
        return 0; /* Singular matrix */
    }
    
    float inv_det = 1.0f / det;
    
    float gp_ho = m12 * m33 - m13 * m32;
    float fp_hn = m11 * m33 - m13 * m31;
    float fo_gn = m11 * m32 - m12 * m31;
    float ep_hm = m10 * m33 - m13 * m30;
    float eo_gm = m10 * m32 - m12 * m30;
    float en_fm = m10 * m31 - m11 * m30;
    
    float gl_hk = m12 * m23 - m13 * m22;
    float fl_hj = m11 * m23 - m13 * m21;
    float fk_gj = m11 * m22 - m12 * m21;
    float el_hi = m10 * m23 - m13 * m20;
    float ek_gi = m10 * m22 - m12 * m20;
    float ej_fi = m10 * m21 - m11 * m20;
    
    out->m[0 * 4 + 0] = a11 * inv_det;
    out->m[0 * 4 + 1] = a12 * inv_det;
    out->m[0 * 4 + 2] = a13 * inv_det;
    out->m[0 * 4 + 3] = a14 * inv_det;
    
    out->m[1 * 4 + 0] = -(m01 * kp_lo - m02 * jp_ln + m03 * jo_kn) * inv_det;
    out->m[1 * 4 + 1] = (m00 * kp_lo - m02 * ip_lm + m03 * io_km) * inv_det;
    out->m[1 * 4 + 2] = -(m00 * jp_ln - m01 * ip_lm + m03 * in_jm) * inv_det;
    out->m[1 * 4 + 3] = (m00 * jo_kn - m01 * io_km + m02 * in_jm) * inv_det;
    
    out->m[2 * 4 + 0] = (m01 * gp_ho - m02 * fp_hn + m03 * fo_gn) * inv_det;
    out->m[2 * 4 + 1] = -(m00 * gp_ho - m02 * ep_hm + m03 * eo_gm) * inv_det;
    out->m[2 * 4 + 2] = (m00 * fp_hn - m01 * ep_hm + m03 * en_fm) * inv_det;
    out->m[2 * 4 + 3] = -(m00 * fo_gn - m01 * eo_gm + m02 * en_fm) * inv_det;
    
    out->m[3 * 4 + 0] = -(m01 * gl_hk - m02 * fl_hj + m03 * fk_gj) * inv_det;
    out->m[3 * 4 + 1] = (m00 * gl_hk - m02 * el_hi + m03 * ek_gi) * inv_det;
    out->m[3 * 4 + 2] = -(m00 * fl_hj - m01 * el_hi + m03 * ej_fi) * inv_det;
    out->m[3 * 4 + 3] = (m00 * fk_gj - m01 * ek_gi + m02 * ej_fi) * inv_det;
    
    return 1;
}

/* ============================================================================
 * Transformation Matrices
 * ============================================================================ */

geomlib_mat4 geomlib_mat4_translate(float x, float y, float z) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    mat.m[3 * 4 + 0] = x;
    mat.m[3 * 4 + 1] = y;
    mat.m[3 * 4 + 2] = z;
    return mat;
}

geomlib_mat4 geomlib_mat4_scale_xyz(float x, float y, float z) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    mat.m[0 * 4 + 0] = x;
    mat.m[1 * 4 + 1] = y;
    mat.m[2 * 4 + 2] = z;
    return mat;
}

geomlib_mat4 geomlib_mat4_rotate_x(float angle_rad) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    float c = cosf(angle_rad);
    float s = sinf(angle_rad);
    mat.m[1 * 4 + 1] = c;
    mat.m[2 * 4 + 1] = -s;
    mat.m[1 * 4 + 2] = s;
    mat.m[2 * 4 + 2] = c;
    return mat;
}

geomlib_mat4 geomlib_mat4_rotate_y(float angle_rad) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    float c = cosf(angle_rad);
    float s = sinf(angle_rad);
    mat.m[0 * 4 + 0] = c;
    mat.m[2 * 4 + 0] = s;
    mat.m[0 * 4 + 2] = -s;
    mat.m[2 * 4 + 2] = c;
    return mat;
}

geomlib_mat4 geomlib_mat4_rotate_z(float angle_rad) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    float c = cosf(angle_rad);
    float s = sinf(angle_rad);
    mat.m[0 * 4 + 0] = c;
    mat.m[1 * 4 + 0] = -s;
    mat.m[0 * 4 + 1] = s;
    mat.m[1 * 4 + 1] = c;
    return mat;
}

geomlib_mat4 geomlib_mat4_rotate_axis(geomlib_vec3 axis, float angle_rad) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    geomlib_vec3 n = geomlib_vec3_normalize(axis);
    float c = cosf(angle_rad);
    float s = sinf(angle_rad);
    float t = 1.0f - c;
    
    mat.m[0 * 4 + 0] = t * n.x * n.x + c;
    mat.m[1 * 4 + 0] = t * n.x * n.y - s * n.z;
    mat.m[2 * 4 + 0] = t * n.x * n.z + s * n.y;
    
    mat.m[0 * 4 + 1] = t * n.x * n.y + s * n.z;
    mat.m[1 * 4 + 1] = t * n.y * n.y + c;
    mat.m[2 * 4 + 1] = t * n.y * n.z - s * n.x;
    
    mat.m[0 * 4 + 2] = t * n.x * n.z - s * n.y;
    mat.m[1 * 4 + 2] = t * n.y * n.z + s * n.x;
    mat.m[2 * 4 + 2] = t * n.z * n.z + c;
    
    return mat;
}

geomlib_mat4 geomlib_mat4_look_at(geomlib_vec3 eye, geomlib_vec3 target, geomlib_vec3 up) {
    geomlib_vec3 f = geomlib_vec3_normalize(geomlib_vec3_sub(target, eye));
    geomlib_vec3 r = geomlib_vec3_normalize(geomlib_vec3_cross(f, up));
    geomlib_vec3 u = geomlib_vec3_cross(r, f);
    
    geomlib_mat4 mat = geomlib_mat4_identity();
    
    mat.m[0 * 4 + 0] = r.x;
    mat.m[1 * 4 + 0] = r.y;
    mat.m[2 * 4 + 0] = r.z;
    mat.m[3 * 4 + 0] = -geomlib_vec3_dot(r, eye);
    
    mat.m[0 * 4 + 1] = u.x;
    mat.m[1 * 4 + 1] = u.y;
    mat.m[2 * 4 + 1] = u.z;
    mat.m[3 * 4 + 1] = -geomlib_vec3_dot(u, eye);
    
    mat.m[0 * 4 + 2] = -f.x;
    mat.m[1 * 4 + 2] = -f.y;
    mat.m[2 * 4 + 2] = -f.z;
    mat.m[3 * 4 + 2] = geomlib_vec3_dot(f, eye);
    
    return mat;
}

geomlib_mat4 geomlib_mat4_perspective(float fov_rad, float aspect, float near_plane, float far_plane) {
    geomlib_mat4 mat = geomlib_mat4_zero();
    float tan_half_fov = tanf(fov_rad / 2.0f);
    
    mat.m[0 * 4 + 0] = 1.0f / (aspect * tan_half_fov);
    mat.m[1 * 4 + 1] = 1.0f / tan_half_fov;
    mat.m[2 * 4 + 2] = -(far_plane + near_plane) / (far_plane - near_plane);
    mat.m[3 * 4 + 2] = -(2.0f * far_plane * near_plane) / (far_plane - near_plane);
    mat.m[2 * 4 + 3] = -1.0f;
    
    return mat;
}

geomlib_mat4 geomlib_mat4_ortho(float left, float right, float bottom, float top, float near_plane, float far_plane) {
    geomlib_mat4 mat = geomlib_mat4_identity();
    
    mat.m[0 * 4 + 0] = 2.0f / (right - left);
    mat.m[1 * 4 + 1] = 2.0f / (top - bottom);
    mat.m[2 * 4 + 2] = -2.0f / (far_plane - near_plane);
    mat.m[3 * 4 + 0] = -(right + left) / (right - left);
    mat.m[3 * 4 + 1] = -(top + bottom) / (top - bottom);
    mat.m[3 * 4 + 2] = -(far_plane + near_plane) / (far_plane - near_plane);
    
    return mat;
}

/* ============================================================================
 * Double-precision Vector2 Operations
 * ============================================================================ */

geomlib_vec2d geomlib_vec2d_create(double x, double y) {
    geomlib_vec2d v;
    v.x = x;
    v.y = y;
    return v;
}

geomlib_vec2d geomlib_vec2d_add(geomlib_vec2d va, geomlib_vec2d vb) {
    geomlib_vec2d v;
    v.x = va.x + vb.x;
    v.y = va.y + vb.y;
    return v;
}

geomlib_vec2d geomlib_vec2d_sub(geomlib_vec2d va, geomlib_vec2d vb) {
    geomlib_vec2d v;
    v.x = va.x - vb.x;
    v.y = va.y - vb.y;
    return v;
}

geomlib_vec2d geomlib_vec2d_mul(geomlib_vec2d v, double s) {
    geomlib_vec2d result;
    result.x = v.x * s;
    result.y = v.y * s;
    return result;
}

double geomlib_vec2d_dot(geomlib_vec2d va, geomlib_vec2d vb) {
    return va.x * vb.x + va.y * vb.y;
}

double geomlib_vec2d_length(geomlib_vec2d v) {
    return sqrt(v.x * v.x + v.y * v.y);
}

geomlib_vec2d geomlib_vec2d_normalize(geomlib_vec2d v) {
    double len = geomlib_vec2d_length(v);
    if (len > 0.0) {
        geomlib_vec2d result;
        double inv_len = 1.0 / len;
        result.x = v.x * inv_len;
        result.y = v.y * inv_len;
        return result;
    }
    return v;
}

/* ============================================================================
 * Double-precision Vector3 Operations
 * ============================================================================ */

geomlib_vec3d geomlib_vec3d_create(double x, double y, double z) {
    geomlib_vec3d v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

geomlib_vec3d geomlib_vec3d_add(geomlib_vec3d va, geomlib_vec3d vb) {
    geomlib_vec3d v;
    v.x = va.x + vb.x;
    v.y = va.y + vb.y;
    v.z = va.z + vb.z;
    return v;
}

geomlib_vec3d geomlib_vec3d_sub(geomlib_vec3d va, geomlib_vec3d vb) {
    geomlib_vec3d v;
    v.x = va.x - vb.x;
    v.y = va.y - vb.y;
    v.z = va.z - vb.z;
    return v;
}

geomlib_vec3d geomlib_vec3d_mul(geomlib_vec3d v, double s) {
    geomlib_vec3d result;
    result.x = v.x * s;
    result.y = v.y * s;
    result.z = v.z * s;
    return result;
}

double geomlib_vec3d_dot(geomlib_vec3d va, geomlib_vec3d vb) {
    return va.x * vb.x + va.y * vb.y + va.z * vb.z;
}

geomlib_vec3d geomlib_vec3d_cross(geomlib_vec3d va, geomlib_vec3d vb) {
    geomlib_vec3d v;
    v.x = va.y * vb.z - va.z * vb.y;
    v.y = va.z * vb.x - va.x * vb.z;
    v.z = va.x * vb.y - va.y * vb.x;
    return v;
}

double geomlib_vec3d_length(geomlib_vec3d v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

geomlib_vec3d geomlib_vec3d_normalize(geomlib_vec3d v) {
    double len = geomlib_vec3d_length(v);
    if (len > 0.0) {
        geomlib_vec3d result;
        double inv_len = 1.0 / len;
        result.x = v.x * inv_len;
        result.y = v.y * inv_len;
        result.z = v.z * inv_len;
        return result;
    }
    return v;
}

/* ============================================================================
 * Double-precision Matrix4 Operations
 * ============================================================================ */

geomlib_mat4d geomlib_mat4d_identity(void) {
    geomlib_mat4d mat;
    geomlib_memset(&mat, 0, sizeof(mat));
    mat.m[0] = 1.0;   /* (0,0) */
    mat.m[5] = 1.0;   /* (1,1) */
    mat.m[10] = 1.0;  /* (2,2) */
    mat.m[15] = 1.0;  /* (3,3) */
    return mat;
}

geomlib_mat4d geomlib_mat4d_mul(geomlib_mat4d ma, geomlib_mat4d mb) {
    geomlib_mat4d mat;
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            mat.m[col * 4 + row] = 
                ma.m[0 * 4 + row] * mb.m[col * 4 + 0] +
                ma.m[1 * 4 + row] * mb.m[col * 4 + 1] +
                ma.m[2 * 4 + row] * mb.m[col * 4 + 2] +
                ma.m[3 * 4 + row] * mb.m[col * 4 + 3];
        }
    }
    return mat;
}

int geomlib_mat4d_inverse(geomlib_mat4d mat, geomlib_mat4d *out) {
    double m00 = mat.m[0 * 4 + 0], m01 = mat.m[1 * 4 + 0], m02 = mat.m[2 * 4 + 0], m03 = mat.m[3 * 4 + 0];
    double m10 = mat.m[0 * 4 + 1], m11 = mat.m[1 * 4 + 1], m12 = mat.m[2 * 4 + 1], m13 = mat.m[3 * 4 + 1];
    double m20 = mat.m[0 * 4 + 2], m21 = mat.m[1 * 4 + 2], m22 = mat.m[2 * 4 + 2], m23 = mat.m[3 * 4 + 2];
    double m30 = mat.m[0 * 4 + 3], m31 = mat.m[1 * 4 + 3], m32 = mat.m[2 * 4 + 3], m33 = mat.m[3 * 4 + 3];
    
    double kp_lo = m22 * m33 - m23 * m32;
    double jp_ln = m21 * m33 - m23 * m31;
    double jo_kn = m21 * m32 - m22 * m31;
    double ip_lm = m20 * m33 - m23 * m30;
    double io_km = m20 * m32 - m22 * m30;
    double in_jm = m20 * m31 - m21 * m30;
    
    double a11 = m11 * kp_lo - m12 * jp_ln + m13 * jo_kn;
    double a12 = -(m10 * kp_lo - m12 * ip_lm + m13 * io_km);
    double a13 = m10 * jp_ln - m11 * ip_lm + m13 * in_jm;
    double a14 = -(m10 * jo_kn - m11 * io_km + m12 * in_jm);
    
    double det = m00 * a11 + m01 * a12 + m02 * a13 + m03 * a14;
    
    if (fabs(det) < 1e-15) {
        return 0; /* Singular matrix */
    }
    
    double inv_det = 1.0 / det;
    
    double gp_ho = m12 * m33 - m13 * m32;
    double fp_hn = m11 * m33 - m13 * m31;
    double fo_gn = m11 * m32 - m12 * m31;
    double ep_hm = m10 * m33 - m13 * m30;
    double eo_gm = m10 * m32 - m12 * m30;
    double en_fm = m10 * m31 - m11 * m30;
    
    double gl_hk = m12 * m23 - m13 * m22;
    double fl_hj = m11 * m23 - m13 * m21;
    double fk_gj = m11 * m22 - m12 * m21;
    double el_hi = m10 * m23 - m13 * m20;
    double ek_gi = m10 * m22 - m12 * m20;
    double ej_fi = m10 * m21 - m11 * m20;
    
    out->m[0 * 4 + 0] = a11 * inv_det;
    out->m[0 * 4 + 1] = a12 * inv_det;
    out->m[0 * 4 + 2] = a13 * inv_det;
    out->m[0 * 4 + 3] = a14 * inv_det;
    
    out->m[1 * 4 + 0] = -(m01 * kp_lo - m02 * jp_ln + m03 * jo_kn) * inv_det;
    out->m[1 * 4 + 1] = (m00 * kp_lo - m02 * ip_lm + m03 * io_km) * inv_det;
    out->m[1 * 4 + 2] = -(m00 * jp_ln - m01 * ip_lm + m03 * in_jm) * inv_det;
    out->m[1 * 4 + 3] = (m00 * jo_kn - m01 * io_km + m02 * in_jm) * inv_det;
    
    out->m[2 * 4 + 0] = (m01 * gp_ho - m02 * fp_hn + m03 * fo_gn) * inv_det;
    out->m[2 * 4 + 1] = -(m00 * gp_ho - m02 * ep_hm + m03 * eo_gm) * inv_det;
    out->m[2 * 4 + 2] = (m00 * fp_hn - m01 * ep_hm + m03 * en_fm) * inv_det;
    out->m[2 * 4 + 3] = -(m00 * fo_gn - m01 * eo_gm + m02 * en_fm) * inv_det;
    
    out->m[3 * 4 + 0] = -(m01 * gl_hk - m02 * fl_hj + m03 * fk_gj) * inv_det;
    out->m[3 * 4 + 1] = (m00 * gl_hk - m02 * el_hi + m03 * ek_gi) * inv_det;
    out->m[3 * 4 + 2] = -(m00 * fl_hj - m01 * el_hi + m03 * ej_fi) * inv_det;
    out->m[3 * 4 + 3] = (m00 * fk_gj - m01 * ek_gi + m02 * ej_fi) * inv_det;
    
    return 1;
}
