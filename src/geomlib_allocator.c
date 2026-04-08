/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Custom allocator implementation.
 */

#include "geomlib/geomlib_allocator.h"
#include "geomlib/geomlib_c11.h"

/* --------------------------------------------------------------------------
 * Default allocator (hosted mode only)
 * -------------------------------------------------------------------------- */

#ifndef GEOMLIB_FREESTANDING

#include <stdlib.h>

static void* default_malloc(size_t size, size_t alignment, void *user_data) {
    (void)user_data;
    (void)alignment; /* stdlib malloc guarantees _Alignof(max_align_t) alignment.
                        Over-aligned requests (alignment > _Alignof(max_align_t))
                        are not supported by this default allocator; provide a
                        custom allocator with aligned_alloc/posix_memalign if needed. */
    return malloc(size);
}

static void* default_realloc(void *ptr, size_t new_size, void *user_data) {
    (void)user_data;
    return realloc(ptr, new_size);
}

static void default_free(void *ptr, void *user_data) {
    (void)user_data;
    free(ptr);
}

static const geomlib_allocator g_default_allocator = {
    default_malloc,
    default_realloc,
    default_free,
    NULL
};

#endif /* !GEOMLIB_FREESTANDING */

/* --------------------------------------------------------------------------
 * Global allocator state
 * -------------------------------------------------------------------------- */

#ifndef GEOMLIB_FREESTANDING
static const geomlib_allocator *g_allocator = &g_default_allocator;
#else
static const geomlib_allocator *g_allocator = NULL;
#endif

void geomlib_set_allocator(const geomlib_allocator *allocator) {
    if (allocator != NULL) {
        g_allocator = allocator;
    } else {
#ifndef GEOMLIB_FREESTANDING
        g_allocator = &g_default_allocator;
#else
        g_allocator = NULL;
#endif
    }
}

const geomlib_allocator* geomlib_get_allocator(void) {
    return g_allocator;
}

#ifndef GEOMLIB_FREESTANDING
const geomlib_allocator* geomlib_default_allocator(void) {
    return &g_default_allocator;
}
#endif

/* --------------------------------------------------------------------------
 * Allocation functions
 * -------------------------------------------------------------------------- */

void* geomlib_malloc(size_t size) {
    if (g_allocator == NULL || g_allocator->malloc_fn == NULL) {
        return NULL;
    }
    return g_allocator->malloc_fn(size, 0, g_allocator->user_data);
}

void* geomlib_aligned_malloc(size_t size, size_t alignment) {
    if (g_allocator == NULL || g_allocator->malloc_fn == NULL) {
        return NULL;
    }
    return g_allocator->malloc_fn(size, alignment, g_allocator->user_data);
}

void* geomlib_calloc(size_t count, size_t size) {
    if (g_allocator == NULL || g_allocator->malloc_fn == NULL) {
        return NULL;
    }

    /* Check for overflow */
    size_t total;
    if (geomlib_safe_mul_size(count, size, &total) != 0) {
        return NULL;
    }

    void *ptr = g_allocator->malloc_fn(total, 0, g_allocator->user_data);
    if (ptr != NULL) {
        geomlib_memset(ptr, 0, total);
    }
    return ptr;
}

void* geomlib_realloc(void *ptr, size_t new_size) {
    if (g_allocator == NULL || g_allocator->realloc_fn == NULL) {
        return NULL;
    }
    return g_allocator->realloc_fn(ptr, new_size, g_allocator->user_data);
}

void geomlib_free(void *ptr) {
    if (ptr == NULL) {
        return;
    }
    if (g_allocator == NULL || g_allocator->free_fn == NULL) {
        return;
    }
    g_allocator->free_fn(ptr, g_allocator->user_data);
}
