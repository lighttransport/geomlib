/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * C11 portability, RAII-like cleanup, and safety utilities.
 *
 * Features:
 * - RAII-like automatic cleanup using __attribute__((cleanup)) (GCC/Clang)
 * - Freestanding memory utilities (memset, memcpy, memmove, memcmp)
 * - Safe arithmetic with overflow detection
 * - Error handling callback
 * - Freestanding math function macros (compiler builtins)
 */

#ifndef GEOMLIB_C11_H
#define GEOMLIB_C11_H

#include "geomlib_allocator.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * Compiler detection
 * ========================================================================== */

#if defined(__GNUC__) || defined(__clang__)
  #define GEOMLIB_HAS_CLEANUP 1
#else
  #define GEOMLIB_HAS_CLEANUP 0
#endif

/* ==========================================================================
 * Freestanding memory utilities
 *
 * In freestanding mode, these provide portable implementations of common
 * memory operations without requiring <string.h>.
 * In hosted mode, they map directly to standard library functions.
 * ========================================================================== */

#ifdef GEOMLIB_FREESTANDING

void* geomlib_memset(void *s, int c, size_t n);
void* geomlib_memcpy(void *dst, const void *src, size_t n);
void* geomlib_memmove(void *dst, const void *src, size_t n);
int   geomlib_memcmp(const void *s1, const void *s2, size_t n);

#else /* Hosted mode: use standard library */

#include <string.h>

#define geomlib_memset  memset
#define geomlib_memcpy  memcpy
#define geomlib_memmove memmove
#define geomlib_memcmp  memcmp

#endif /* GEOMLIB_FREESTANDING */

/* ==========================================================================
 * Freestanding stdlib utilities (qsort, abs)
 *
 * In freestanding mode, these provide portable implementations.
 * In hosted mode, they map to standard library functions.
 * ========================================================================== */

#ifdef GEOMLIB_FREESTANDING

/** Portable qsort implementation (Shell sort) for freestanding environments. */
void geomlib_qsort(void *base, size_t nmemb, size_t size,
                   int (*compar)(const void *, const void *));

/** Portable abs for freestanding environments. */
static inline int geomlib_abs(int x) { return x < 0 ? -x : x; }

#else /* Hosted mode */

#include <stdlib.h>
#define geomlib_qsort  qsort
#define geomlib_abs    abs

#endif /* GEOMLIB_FREESTANDING */

/* ==========================================================================
 * Freestanding math macros
 *
 * In freestanding mode, math functions are mapped to compiler builtins
 * (supported by GCC and Clang). Users can override these by defining
 * the macros before including this header.
 *
 * In hosted mode, standard <math.h> is used.
 * ========================================================================== */

#ifdef GEOMLIB_FREESTANDING

/* Float versions */
#ifndef geomlib_sqrtf
  #define geomlib_sqrtf(x)  __builtin_sqrtf(x)
#endif
#ifndef geomlib_fabsf
  #define geomlib_fabsf(x)  __builtin_fabsf(x)
#endif
#ifndef geomlib_sinf
  #define geomlib_sinf(x)   __builtin_sinf(x)
#endif
#ifndef geomlib_cosf
  #define geomlib_cosf(x)   __builtin_cosf(x)
#endif
#ifndef geomlib_tanf
  #define geomlib_tanf(x)   __builtin_tanf(x)
#endif
#ifndef geomlib_asinf
  #define geomlib_asinf(x)  __builtin_asinf(x)
#endif
#ifndef geomlib_acosf
  #define geomlib_acosf(x)  __builtin_acosf(x)
#endif
#ifndef geomlib_atan2f
  #define geomlib_atan2f(y,x) __builtin_atan2f(y,x)
#endif
#ifndef geomlib_floorf
  #define geomlib_floorf(x) __builtin_floorf(x)
#endif
#ifndef geomlib_ceilf
  #define geomlib_ceilf(x)  __builtin_ceilf(x)
#endif
#ifndef geomlib_fmodf
  #define geomlib_fmodf(x,y) __builtin_fmodf(x,y)
#endif
#ifndef geomlib_powf
  #define geomlib_powf(x,y) __builtin_powf(x,y)
#endif
#ifndef geomlib_expf
  #define geomlib_expf(x)   __builtin_expf(x)
#endif
#ifndef geomlib_logf
  #define geomlib_logf(x)   __builtin_logf(x)
#endif
#ifndef geomlib_log2f
  #define geomlib_log2f(x)  __builtin_log2f(x)
#endif
#ifndef geomlib_cbrtf
  #define geomlib_cbrtf(x)  __builtin_cbrtf(x)
#endif
#ifndef geomlib_fmaxf
  #define geomlib_fmaxf(x,y) __builtin_fmaxf(x,y)
#endif
#ifndef geomlib_fminf
  #define geomlib_fminf(x,y) __builtin_fminf(x,y)
#endif
#ifndef geomlib_roundf
  #define geomlib_roundf(x) __builtin_roundf(x)
#endif

/* Double versions */
#ifndef geomlib_sqrt
  #define geomlib_sqrt(x)   __builtin_sqrt(x)
#endif
#ifndef geomlib_fabs
  #define geomlib_fabs(x)   __builtin_fabs(x)
#endif
#ifndef geomlib_sin
  #define geomlib_sin(x)    __builtin_sin(x)
#endif
#ifndef geomlib_cos
  #define geomlib_cos(x)    __builtin_cos(x)
#endif
#ifndef geomlib_tan
  #define geomlib_tan(x)    __builtin_tan(x)
#endif
#ifndef geomlib_asin
  #define geomlib_asin(x)   __builtin_asin(x)
#endif
#ifndef geomlib_acos
  #define geomlib_acos(x)   __builtin_acos(x)
#endif
#ifndef geomlib_atan2
  #define geomlib_atan2(y,x) __builtin_atan2(y,x)
#endif
#ifndef geomlib_floor
  #define geomlib_floor(x)  __builtin_floor(x)
#endif
#ifndef geomlib_ceil
  #define geomlib_ceil(x)   __builtin_ceil(x)
#endif
#ifndef geomlib_fmod
  #define geomlib_fmod(x,y) __builtin_fmod(x,y)
#endif
#ifndef geomlib_pow
  #define geomlib_pow(x,y)  __builtin_pow(x,y)
#endif
#ifndef geomlib_exp
  #define geomlib_exp(x)    __builtin_exp(x)
#endif
#ifndef geomlib_log
  #define geomlib_log(x)    __builtin_log(x)
#endif
#ifndef geomlib_log2
  #define geomlib_log2(x)   __builtin_log2(x)
#endif
#ifndef geomlib_cbrt
  #define geomlib_cbrt(x)   __builtin_cbrt(x)
#endif
#ifndef geomlib_fmax
  #define geomlib_fmax(x,y) __builtin_fmax(x,y)
#endif
#ifndef geomlib_fmin
  #define geomlib_fmin(x,y) __builtin_fmin(x,y)
#endif
#ifndef geomlib_round
  #define geomlib_round(x)  __builtin_round(x)
#endif

/* Constants */
#ifndef GEOMLIB_PI
  #define GEOMLIB_PI  3.14159265358979323846
#endif
#ifndef GEOMLIB_PIf
  #define GEOMLIB_PIf 3.14159265358979323846f
#endif

#else /* Hosted mode: use standard math */

#include <math.h>

#define geomlib_sqrtf  sqrtf
#define geomlib_fabsf  fabsf
#define geomlib_sinf   sinf
#define geomlib_cosf   cosf
#define geomlib_tanf   tanf
#define geomlib_asinf  asinf
#define geomlib_acosf  acosf
#define geomlib_atan2f atan2f
#define geomlib_floorf floorf
#define geomlib_ceilf  ceilf
#define geomlib_fmodf  fmodf
#define geomlib_powf   powf
#define geomlib_expf   expf
#define geomlib_logf   logf
#define geomlib_log2f  log2f
#define geomlib_cbrtf  cbrtf
#define geomlib_fmaxf  fmaxf
#define geomlib_fminf  fminf
#define geomlib_roundf roundf

#define geomlib_sqrt   sqrt
#define geomlib_fabs   fabs
#define geomlib_sin    sin
#define geomlib_cos    cos
#define geomlib_tan    tan
#define geomlib_asin   asin
#define geomlib_acos   acos
#define geomlib_atan2  atan2
#define geomlib_floor  floor
#define geomlib_ceil   ceil
#define geomlib_fmod   fmod
#define geomlib_pow    pow
#define geomlib_exp    exp
#define geomlib_log    log
#define geomlib_log2   log2
#define geomlib_cbrt   cbrt
#define geomlib_fmax   fmax
#define geomlib_fmin   fmin
#define geomlib_round  round

#ifndef GEOMLIB_PI
  #define GEOMLIB_PI   3.14159265358979323846
#endif
#ifndef GEOMLIB_PIf
  #define GEOMLIB_PIf  3.14159265358979323846f
#endif

#endif /* GEOMLIB_FREESTANDING */

/* ==========================================================================
 * RAII-like automatic cleanup (GCC/Clang only)
 *
 * Usage:
 *   {
 *       GEOMLIB_AUTO_TRIMESH geomlib_trimesh mesh;
 *       geomlib_trimesh_init(&mesh, 100, 100);
 *       // ... use mesh ...
 *   } // mesh is automatically freed here
 *
 * For heap-allocated (pointer) types:
 *   {
 *       GEOMLIB_AUTO_FREE void *ptr = geomlib_malloc(1024);
 *       // ... use ptr ...
 *   } // ptr is automatically freed here
 * ========================================================================== */

/* Forward declarations of library types */
struct geomlib_trimesh;
struct geomlib_polymesh;
struct geomlib_polygon2d;
struct geomlib_polygon3d;
struct geomlib_bvh;
struct geomlib_mesh_sample_distribution;

/* Cleanup functions - free the resource pointed to by the stack variable */
void geomlib_auto_free_trimesh(struct geomlib_trimesh *p);
void geomlib_auto_free_polymesh(struct geomlib_polymesh *p);
void geomlib_auto_free_polygon2d(struct geomlib_polygon2d *p);
void geomlib_auto_free_polygon3d(struct geomlib_polygon3d *p);
void geomlib_auto_free_bvh(struct geomlib_bvh *p);
void geomlib_auto_free_sample_dist(struct geomlib_mesh_sample_distribution *p);

/* Generic pointer cleanup - frees via geomlib_free */
void geomlib_auto_free_ptr(void **p);

#if GEOMLIB_HAS_CLEANUP

  #define GEOMLIB_AUTO_TRIMESH     __attribute__((cleanup(geomlib_auto_free_trimesh)))
  #define GEOMLIB_AUTO_POLYMESH    __attribute__((cleanup(geomlib_auto_free_polymesh)))
  #define GEOMLIB_AUTO_POLYGON2D   __attribute__((cleanup(geomlib_auto_free_polygon2d)))
  #define GEOMLIB_AUTO_POLYGON3D   __attribute__((cleanup(geomlib_auto_free_polygon3d)))
  #define GEOMLIB_AUTO_BVH         __attribute__((cleanup(geomlib_auto_free_bvh)))
  #define GEOMLIB_AUTO_SAMPLE_DIST __attribute__((cleanup(geomlib_auto_free_sample_dist)))
  #define GEOMLIB_AUTO_FREE        __attribute__((cleanup(geomlib_auto_free_ptr)))

#else /* Compiler without cleanup attribute support */

  #define GEOMLIB_AUTO_TRIMESH
  #define GEOMLIB_AUTO_POLYMESH
  #define GEOMLIB_AUTO_POLYGON2D
  #define GEOMLIB_AUTO_POLYGON3D
  #define GEOMLIB_AUTO_BVH
  #define GEOMLIB_AUTO_SAMPLE_DIST
  #define GEOMLIB_AUTO_FREE

#endif /* GEOMLIB_HAS_CLEANUP */

/* ==========================================================================
 * Safe arithmetic (overflow detection)
 * ========================================================================== */

/**
 * Safe multiply for size_t. Returns 0 on success, non-zero on overflow.
 * Result is stored in *result.
 */
int geomlib_safe_mul_size(size_t a, size_t b, size_t *result);

/**
 * Safe add for size_t. Returns 0 on success, non-zero on overflow.
 * Result is stored in *result.
 */
int geomlib_safe_add_size(size_t a, size_t b, size_t *result);

/* ==========================================================================
 * Error handling
 * ========================================================================== */

/** Error severity levels */
typedef enum geomlib_error_level {
    GEOMLIB_ERROR_WARNING = 0,
    GEOMLIB_ERROR_ERROR   = 1,
    GEOMLIB_ERROR_FATAL   = 2
} geomlib_error_level;

/**
 * Error handler callback type.
 *
 * @param level     Severity of the error
 * @param message   Human-readable error description
 * @param file      Source file where the error occurred (may be NULL)
 * @param line      Source line where the error occurred
 * @param user_data Opaque user pointer set with geomlib_set_error_handler
 */
typedef void (*geomlib_error_handler_fn)(
    geomlib_error_level level,
    const char *message,
    const char *file,
    int line,
    void *user_data
);

/**
 * Set a custom error handler.
 * Pass NULL to disable error reporting (default).
 */
void geomlib_set_error_handler(geomlib_error_handler_fn handler, void *user_data);

/**
 * Report an error through the registered handler.
 * No-op if no handler is set.
 */
void geomlib_report_error(geomlib_error_level level, const char *message,
                          const char *file, int line);

/** Convenience macro for reporting errors with automatic file/line */
#define GEOMLIB_ERROR(level, msg) \
    geomlib_report_error((level), (msg), __FILE__, __LINE__)

#define GEOMLIB_WARN(msg)  GEOMLIB_ERROR(GEOMLIB_ERROR_WARNING, (msg))
#define GEOMLIB_ERR(msg)   GEOMLIB_ERROR(GEOMLIB_ERROR_ERROR, (msg))
#define GEOMLIB_FATAL(msg) GEOMLIB_ERROR(GEOMLIB_ERROR_FATAL, (msg))

/* ==========================================================================
 * Null-safe pointer validation
 * ========================================================================== */

/** Check if pointer is non-NULL. Returns the pointer or NULL. */
#define GEOMLIB_CHECK_PTR(ptr) ((ptr) != NULL ? (ptr) : NULL)

/** Return error code if pointer is NULL */
#define GEOMLIB_REQUIRE_NON_NULL(ptr, retval) \
    do { \
        if ((ptr) == NULL) { \
            GEOMLIB_ERR("NULL pointer argument"); \
            return (retval); \
        } \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_C11_H */
