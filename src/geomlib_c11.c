/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * C11 utilities: RAII cleanup, freestanding memory ops, safe arithmetic,
 * error handling.
 */

#include "geomlib/geomlib_c11.h"
#include "geomlib/geomlib_allocator.h"

/* Include full type definitions for cleanup functions */
#include "geomlib/geomlib_mesh.h"
#include "geomlib/geomlib_bvh.h"
#include "geomlib/geomlib_montecarlo.h"

/* ==========================================================================
 * Freestanding memory utilities
 * ========================================================================== */

#ifdef GEOMLIB_FREESTANDING

void* geomlib_memset(void *s, int c, size_t n) {
    unsigned char *p = (unsigned char *)s;
    unsigned char val = (unsigned char)c;
    while (n--) {
        *p++ = val;
    }
    return s;
}

void* geomlib_memcpy(void *dst, const void *src, size_t n) {
    unsigned char *d = (unsigned char *)dst;
    const unsigned char *s = (const unsigned char *)src;
    while (n--) {
        *d++ = *s++;
    }
    return dst;
}

void* geomlib_memmove(void *dst, const void *src, size_t n) {
    unsigned char *d = (unsigned char *)dst;
    const unsigned char *s = (const unsigned char *)src;
    if (d < s) {
        while (n--) {
            *d++ = *s++;
        }
    } else if (d > s) {
        d += n;
        s += n;
        while (n--) {
            *--d = *--s;
        }
    }
    return dst;
}

int geomlib_memcmp(const void *s1, const void *s2, size_t n) {
    const unsigned char *p1 = (const unsigned char *)s1;
    const unsigned char *p2 = (const unsigned char *)s2;
    while (n--) {
        if (*p1 != *p2) {
            return (int)*p1 - (int)*p2;
        }
        p1++;
        p2++;
    }
    return 0;
}

#endif /* GEOMLIB_FREESTANDING */

/* ==========================================================================
 * RAII cleanup functions
 * ========================================================================== */

void geomlib_auto_free_trimesh(struct geomlib_trimesh *p) {
    if (p != NULL) {
        geomlib_trimesh_free(p);
    }
}

void geomlib_auto_free_polymesh(struct geomlib_polymesh *p) {
    if (p != NULL) {
        geomlib_polymesh_free(p);
    }
}

void geomlib_auto_free_polygon2d(struct geomlib_polygon2d *p) {
    if (p != NULL) {
        geomlib_polygon2d_free(p);
    }
}

void geomlib_auto_free_polygon3d(struct geomlib_polygon3d *p) {
    if (p != NULL) {
        geomlib_polygon3d_free(p);
    }
}

void geomlib_auto_free_bvh(struct geomlib_bvh *p) {
    if (p != NULL) {
        geomlib_bvh_free(p);
    }
}

void geomlib_auto_free_sample_dist(struct geomlib_mesh_sample_distribution *p) {
    if (p != NULL) {
        geomlib_mesh_free_sample_distribution(p);
    }
}

void geomlib_auto_free_ptr(void **p) {
    if (p != NULL && *p != NULL) {
        geomlib_free(*p);
        *p = NULL;
    }
}

/* ==========================================================================
 * Safe arithmetic
 * ========================================================================== */

int geomlib_safe_mul_size(size_t a, size_t b, size_t *result) {
    if (result == NULL) {
        return -1;
    }
    if (a == 0 || b == 0) {
        *result = 0;
        return 0;
    }
    /* Check for overflow: a * b > SIZE_MAX iff a > SIZE_MAX / b */
    if (a > (size_t)-1 / b) {
        *result = 0;
        return -1;
    }
    *result = a * b;
    return 0;
}

int geomlib_safe_add_size(size_t a, size_t b, size_t *result) {
    if (result == NULL) {
        return -1;
    }
    if (a > (size_t)-1 - b) {
        *result = 0;
        return -1;
    }
    *result = a + b;
    return 0;
}

/* ==========================================================================
 * Error handling
 * ========================================================================== */

static geomlib_error_handler_fn g_error_handler = NULL;
static void *g_error_user_data = NULL;

void geomlib_set_error_handler(geomlib_error_handler_fn handler, void *user_data) {
    g_error_handler = handler;
    g_error_user_data = user_data;
}

void geomlib_report_error(geomlib_error_level level, const char *message,
                          const char *file, int line) {
    if (g_error_handler != NULL) {
        g_error_handler(level, message, file, line, g_error_user_data);
    }
}
