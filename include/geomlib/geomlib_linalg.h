/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Linear Algebra Foundation
 * Core functions implemented in C for maximum portability
 */

#ifndef GEOMLIB_LINALG_H
#define GEOMLIB_LINALG_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Vector Types
 * ============================================================================ */

typedef struct geomlib_vec2 {
    float x, y;
} geomlib_vec2;

typedef struct geomlib_vec3 {
    float x, y, z;
} geomlib_vec3;

typedef struct geomlib_vec4 {
    float x, y, z, w;
} geomlib_vec4;

typedef struct geomlib_vec2d {
    double x, y;
} geomlib_vec2d;

typedef struct geomlib_vec3d {
    double x, y, z;
} geomlib_vec3d;

typedef struct geomlib_vec4d {
    double x, y, z, w;
} geomlib_vec4d;

/* ============================================================================
 * Matrix Types (Column-major order)
 * ============================================================================ */

typedef struct geomlib_mat3 {
    float m[9];  /* Column-major: m[col*3 + row] */
} geomlib_mat3;

typedef struct geomlib_mat4 {
    float m[16]; /* Column-major: m[col*4 + row] */
} geomlib_mat4;

typedef struct geomlib_mat3d {
    double m[9];
} geomlib_mat3d;

typedef struct geomlib_mat4d {
    double m[16];
} geomlib_mat4d;

/* ============================================================================
 * Vector2 Operations (float)
 * ============================================================================ */

geomlib_vec2 geomlib_vec2_create(float x, float y);
geomlib_vec2 geomlib_vec2_add(geomlib_vec2 a, geomlib_vec2 b);
geomlib_vec2 geomlib_vec2_sub(geomlib_vec2 a, geomlib_vec2 b);
geomlib_vec2 geomlib_vec2_mul(geomlib_vec2 v, float s);
geomlib_vec2 geomlib_vec2_div(geomlib_vec2 v, float s);
float geomlib_vec2_dot(geomlib_vec2 a, geomlib_vec2 b);
float geomlib_vec2_length(geomlib_vec2 v);
float geomlib_vec2_length_squared(geomlib_vec2 v);
geomlib_vec2 geomlib_vec2_normalize(geomlib_vec2 v);
geomlib_vec2 geomlib_vec2_negate(geomlib_vec2 v);

/* ============================================================================
 * Vector3 Operations (float)
 * ============================================================================ */

geomlib_vec3 geomlib_vec3_create(float x, float y, float z);
geomlib_vec3 geomlib_vec3_add(geomlib_vec3 a, geomlib_vec3 b);
geomlib_vec3 geomlib_vec3_sub(geomlib_vec3 a, geomlib_vec3 b);
geomlib_vec3 geomlib_vec3_mul(geomlib_vec3 v, float s);
geomlib_vec3 geomlib_vec3_div(geomlib_vec3 v, float s);
float geomlib_vec3_dot(geomlib_vec3 a, geomlib_vec3 b);
geomlib_vec3 geomlib_vec3_cross(geomlib_vec3 a, geomlib_vec3 b);
float geomlib_vec3_length(geomlib_vec3 v);
float geomlib_vec3_length_squared(geomlib_vec3 v);
geomlib_vec3 geomlib_vec3_normalize(geomlib_vec3 v);
geomlib_vec3 geomlib_vec3_negate(geomlib_vec3 v);
geomlib_vec3 geomlib_vec3_lerp(geomlib_vec3 a, geomlib_vec3 b, float t);

/* ============================================================================
 * Vector4 Operations (float)
 * ============================================================================ */

geomlib_vec4 geomlib_vec4_create(float x, float y, float z, float w);
geomlib_vec4 geomlib_vec4_add(geomlib_vec4 a, geomlib_vec4 b);
geomlib_vec4 geomlib_vec4_sub(geomlib_vec4 a, geomlib_vec4 b);
geomlib_vec4 geomlib_vec4_mul(geomlib_vec4 v, float s);
geomlib_vec4 geomlib_vec4_div(geomlib_vec4 v, float s);
float geomlib_vec4_dot(geomlib_vec4 a, geomlib_vec4 b);
float geomlib_vec4_length(geomlib_vec4 v);
float geomlib_vec4_length_squared(geomlib_vec4 v);
geomlib_vec4 geomlib_vec4_normalize(geomlib_vec4 v);
geomlib_vec4 geomlib_vec4_negate(geomlib_vec4 v);

/* ============================================================================
 * Matrix3 Operations (float)
 * ============================================================================ */

geomlib_mat3 geomlib_mat3_identity(void);
geomlib_mat3 geomlib_mat3_zero(void);
geomlib_mat3 geomlib_mat3_add(geomlib_mat3 a, geomlib_mat3 b);
geomlib_mat3 geomlib_mat3_sub(geomlib_mat3 a, geomlib_mat3 b);
geomlib_mat3 geomlib_mat3_mul(geomlib_mat3 a, geomlib_mat3 b);
geomlib_mat3 geomlib_mat3_scale(geomlib_mat3 m, float s);
geomlib_vec3 geomlib_mat3_mul_vec3(geomlib_mat3 m, geomlib_vec3 v);
geomlib_mat3 geomlib_mat3_transpose(geomlib_mat3 m);
float geomlib_mat3_determinant(geomlib_mat3 m);
int geomlib_mat3_inverse(geomlib_mat3 m, geomlib_mat3 *out);

/* ============================================================================
 * Matrix4 Operations (float)
 * ============================================================================ */

geomlib_mat4 geomlib_mat4_identity(void);
geomlib_mat4 geomlib_mat4_zero(void);
geomlib_mat4 geomlib_mat4_add(geomlib_mat4 a, geomlib_mat4 b);
geomlib_mat4 geomlib_mat4_sub(geomlib_mat4 a, geomlib_mat4 b);
geomlib_mat4 geomlib_mat4_mul(geomlib_mat4 a, geomlib_mat4 b);
geomlib_mat4 geomlib_mat4_scale(geomlib_mat4 m, float s);
geomlib_vec4 geomlib_mat4_mul_vec4(geomlib_mat4 m, geomlib_vec4 v);
geomlib_mat4 geomlib_mat4_transpose(geomlib_mat4 m);
float geomlib_mat4_determinant(geomlib_mat4 m);
int geomlib_mat4_inverse(geomlib_mat4 m, geomlib_mat4 *out);

/* Transformation matrices */
geomlib_mat4 geomlib_mat4_translate(float x, float y, float z);
geomlib_mat4 geomlib_mat4_scale_xyz(float x, float y, float z);
geomlib_mat4 geomlib_mat4_rotate_x(float angle_rad);
geomlib_mat4 geomlib_mat4_rotate_y(float angle_rad);
geomlib_mat4 geomlib_mat4_rotate_z(float angle_rad);
geomlib_mat4 geomlib_mat4_rotate_axis(geomlib_vec3 axis, float angle_rad);
geomlib_mat4 geomlib_mat4_look_at(geomlib_vec3 eye, geomlib_vec3 target, geomlib_vec3 up);
geomlib_mat4 geomlib_mat4_perspective(float fov_rad, float aspect, float near_plane, float far_plane);
geomlib_mat4 geomlib_mat4_ortho(float left, float right, float bottom, float top, float near_plane, float far_plane);

/* ============================================================================
 * Double-precision variants (selected operations)
 * ============================================================================ */

geomlib_vec3d geomlib_vec3d_create(double x, double y, double z);
geomlib_vec3d geomlib_vec3d_add(geomlib_vec3d a, geomlib_vec3d b);
geomlib_vec3d geomlib_vec3d_sub(geomlib_vec3d a, geomlib_vec3d b);
geomlib_vec3d geomlib_vec3d_mul(geomlib_vec3d v, double s);
double geomlib_vec3d_dot(geomlib_vec3d a, geomlib_vec3d b);
geomlib_vec3d geomlib_vec3d_cross(geomlib_vec3d a, geomlib_vec3d b);
double geomlib_vec3d_length(geomlib_vec3d v);
geomlib_vec3d geomlib_vec3d_normalize(geomlib_vec3d v);

geomlib_mat4d geomlib_mat4d_identity(void);
geomlib_mat4d geomlib_mat4d_mul(geomlib_mat4d a, geomlib_mat4d b);
int geomlib_mat4d_inverse(geomlib_mat4d m, geomlib_mat4d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_LINALG_H */
