/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * BVH and Ray Tracing Implementation
 */

#include "geomlib/geomlib_bvh.h"
#include "geomlib/geomlib_montecarlo.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define GEOMLIB_BVH_MAX_PRIMS_IN_LEAF 4
#define GEOMLIB_BVH_SAH_BUCKETS 12
#define GEOMLIB_BVH_TRAVERSAL_COST 0.125f
#define GEOMLIB_BVH_INTERSECTION_COST 1.0f

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

geomlib_ray geomlib_ray_create(geomlib_vec3 origin, geomlib_vec3 direction) {
    geomlib_ray ray;
    ray.origin = origin;
    ray.direction = geomlib_vec3_normalize(direction);
    ray.t_min = 0.0f;
    ray.t_max = FLT_MAX;
    return ray;
}

geomlib_vec3 geomlib_ray_at(const geomlib_ray *ray, float t) {
    geomlib_vec3 result;
    result.x = ray->origin.x + t * ray->direction.x;
    result.y = ray->origin.y + t * ray->direction.y;
    result.z = ray->origin.z + t * ray->direction.z;
    return result;
}

/* ============================================================================
 * Ray-Primitive Intersection
 * ============================================================================ */

bool geomlib_ray_triangle_intersect(
    const geomlib_ray *ray,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float *t_out,
    float *u_out,
    float *v_out
) {
    /* Möller-Trumbore ray-triangle intersection algorithm */
    const float EPSILON = 1e-8f;

    geomlib_vec3 edge1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 edge2 = geomlib_vec3_sub(v2, v0);

    geomlib_vec3 h = geomlib_vec3_cross(ray->direction, edge2);
    float a = geomlib_vec3_dot(edge1, h);

    /* Ray parallel to triangle */
    if (fabsf(a) < EPSILON) {
        return false;
    }

    float f = 1.0f / a;
    geomlib_vec3 s = geomlib_vec3_sub(ray->origin, v0);
    float u = f * geomlib_vec3_dot(s, h);

    if (u < 0.0f || u > 1.0f) {
        return false;
    }

    geomlib_vec3 q = geomlib_vec3_cross(s, edge1);
    float v = f * geomlib_vec3_dot(ray->direction, q);

    if (v < 0.0f || u + v > 1.0f) {
        return false;
    }

    float t = f * geomlib_vec3_dot(edge2, q);

    if (t < ray->t_min || t > ray->t_max) {
        return false;
    }

    if (t_out) *t_out = t;
    if (u_out) *u_out = u;
    if (v_out) *v_out = v;

    return true;
}

bool geomlib_ray_aabb_intersect(
    const geomlib_ray *ray,
    const geomlib_aabb *aabb
) {
    /* Optimized slab method */
    float inv_dir_x = 1.0f / ray->direction.x;
    float inv_dir_y = 1.0f / ray->direction.y;
    float inv_dir_z = 1.0f / ray->direction.z;

    float t1 = (aabb->min.x - ray->origin.x) * inv_dir_x;
    float t2 = (aabb->max.x - ray->origin.x) * inv_dir_x;
    float t3 = (aabb->min.y - ray->origin.y) * inv_dir_y;
    float t4 = (aabb->max.y - ray->origin.y) * inv_dir_y;
    float t5 = (aabb->min.z - ray->origin.z) * inv_dir_z;
    float t6 = (aabb->max.z - ray->origin.z) * inv_dir_z;

    float tmin = fmaxf(fmaxf(fminf(t1, t2), fminf(t3, t4)), fminf(t5, t6));
    float tmax = fminf(fminf(fmaxf(t1, t2), fmaxf(t3, t4)), fmaxf(t5, t6));

    /* Ray intersects AABB if tmax >= max(0, tmin) && tmin < ray.t_max */
    return tmax >= fmaxf(0.0f, tmin) && tmin < ray->t_max;
}

/* ============================================================================
 * BVH Construction Helpers
 * ============================================================================ */

typedef struct bvh_build_node {
    geomlib_aabb bounds;
    uint32_t first_prim_offset;
    uint32_t prim_count;
    uint32_t split_axis;
    uint32_t left_child;
    uint32_t right_child;
} bvh_build_node;

typedef struct bvh_primitive_info {
    uint32_t primitive_index;
    geomlib_vec3 centroid;
    geomlib_aabb bounds;
} bvh_primitive_info;

typedef struct sah_bucket {
    uint32_t count;
    geomlib_aabb bounds;
} sah_bucket;

static geomlib_aabb compute_triangle_bounds(
    const geomlib_trimesh *mesh,
    uint32_t tri_idx
) {
    uint32_t i0 = mesh->indices[tri_idx * 3 + 0];
    uint32_t i1 = mesh->indices[tri_idx * 3 + 1];
    uint32_t i2 = mesh->indices[tri_idx * 3 + 2];

    geomlib_vec3 v0 = mesh->vertices[i0];
    geomlib_vec3 v1 = mesh->vertices[i1];
    geomlib_vec3 v2 = mesh->vertices[i2];

    geomlib_aabb bounds;
    bounds.min.x = fminf(fminf(v0.x, v1.x), v2.x);
    bounds.min.y = fminf(fminf(v0.y, v1.y), v2.y);
    bounds.min.z = fminf(fminf(v0.z, v1.z), v2.z);
    bounds.max.x = fmaxf(fmaxf(v0.x, v1.x), v2.x);
    bounds.max.y = fmaxf(fmaxf(v0.y, v1.y), v2.y);
    bounds.max.z = fmaxf(fmaxf(v0.z, v1.z), v2.z);

    return bounds;
}

static geomlib_vec3 aabb_centroid(const geomlib_aabb *aabb) {
    geomlib_vec3 centroid;
    centroid.x = (aabb->min.x + aabb->max.x) * 0.5f;
    centroid.y = (aabb->min.y + aabb->max.y) * 0.5f;
    centroid.z = (aabb->min.z + aabb->max.z) * 0.5f;
    return centroid;
}

static float aabb_surface_area(const geomlib_aabb *aabb) {
    geomlib_vec3 size = geomlib_vec3_sub(aabb->max, aabb->min);
    return 2.0f * (size.x * size.y + size.y * size.z + size.z * size.x);
}

static geomlib_aabb aabb_union(const geomlib_aabb *a, const geomlib_aabb *b) {
    geomlib_aabb result;
    result.min.x = fminf(a->min.x, b->min.x);
    result.min.y = fminf(a->min.y, b->min.y);
    result.min.z = fminf(a->min.z, b->min.z);
    result.max.x = fmaxf(a->max.x, b->max.x);
    result.max.y = fmaxf(a->max.y, b->max.y);
    result.max.z = fmaxf(a->max.z, b->max.z);
    return result;
}

static int compare_primitive_x(const void *a, const void *b) {
    const bvh_primitive_info *pa = (const bvh_primitive_info *)a;
    const bvh_primitive_info *pb = (const bvh_primitive_info *)b;
    if (pa->centroid.x < pb->centroid.x) return -1;
    if (pa->centroid.x > pb->centroid.x) return 1;
    return 0;
}

static int compare_primitive_y(const void *a, const void *b) {
    const bvh_primitive_info *pa = (const bvh_primitive_info *)a;
    const bvh_primitive_info *pb = (const bvh_primitive_info *)b;
    if (pa->centroid.y < pb->centroid.y) return -1;
    if (pa->centroid.y > pb->centroid.y) return 1;
    return 0;
}

static int compare_primitive_z(const void *a, const void *b) {
    const bvh_primitive_info *pa = (const bvh_primitive_info *)a;
    const bvh_primitive_info *pb = (const bvh_primitive_info *)b;
    if (pa->centroid.z < pb->centroid.z) return -1;
    if (pa->centroid.z > pb->centroid.z) return 1;
    return 0;
}

static uint32_t build_bvh_recursive(
    geomlib_bvh *bvh,
    bvh_primitive_info *primitives,
    uint32_t start,
    uint32_t end,
    geomlib_bvh_build_quality quality
) {
    uint32_t node_idx = bvh->node_count++;

    if (bvh->node_count > bvh->node_capacity) {
        /* Expand node array */
        bvh->node_capacity *= 2;
        bvh->nodes = (geomlib_bvh_node *)realloc(
            bvh->nodes,
            sizeof(geomlib_bvh_node) * bvh->node_capacity
        );
    }

    geomlib_bvh_node *node = &bvh->nodes[node_idx];
    uint32_t prim_count = end - start;

    /* Compute bounds for all primitives */
    geomlib_aabb bounds = primitives[start].bounds;
    geomlib_aabb centroid_bounds = {{primitives[start].centroid.x, primitives[start].centroid.y, primitives[start].centroid.z},
                                     {primitives[start].centroid.x, primitives[start].centroid.y, primitives[start].centroid.z}};

    for (uint32_t i = start + 1; i < end; i++) {
        bounds = aabb_union(&bounds, &primitives[i].bounds);
        centroid_bounds.min.x = fminf(centroid_bounds.min.x, primitives[i].centroid.x);
        centroid_bounds.min.y = fminf(centroid_bounds.min.y, primitives[i].centroid.y);
        centroid_bounds.min.z = fminf(centroid_bounds.min.z, primitives[i].centroid.z);
        centroid_bounds.max.x = fmaxf(centroid_bounds.max.x, primitives[i].centroid.x);
        centroid_bounds.max.y = fmaxf(centroid_bounds.max.y, primitives[i].centroid.y);
        centroid_bounds.max.z = fmaxf(centroid_bounds.max.z, primitives[i].centroid.z);
    }

    node->bounds = bounds;

    /* Create leaf if few primitives or centroids degenerate */
    geomlib_vec3 centroid_extent = geomlib_vec3_sub(centroid_bounds.max, centroid_bounds.min);
    bool degenerate = (centroid_extent.x < 1e-6f && centroid_extent.y < 1e-6f && centroid_extent.z < 1e-6f);

    if (prim_count <= GEOMLIB_BVH_MAX_PRIMS_IN_LEAF || degenerate) {
        /* Create leaf node */
        node->is_leaf = true;
        node->data.leaf.primitive_offset = start;
        node->data.leaf.primitive_count = prim_count;
        node->axis = 0;
        return node_idx;
    }

    /* Choose split axis (longest centroid extent) */
    uint32_t split_axis = 0;
    if (centroid_extent.y > centroid_extent.x) split_axis = 1;
    if (centroid_extent.z > centroid_extent.x && centroid_extent.z > centroid_extent.y) split_axis = 2;

    node->axis = split_axis;
    node->is_leaf = false;

    uint32_t mid = start + prim_count / 2;

    /* Split using SAH or median based on quality */
    if (quality >= GEOMLIB_BVH_QUALITY_MEDIUM && prim_count > GEOMLIB_BVH_MAX_PRIMS_IN_LEAF) {
        /* SAH split */
        sah_bucket buckets[GEOMLIB_BVH_SAH_BUCKETS];
        memset(buckets, 0, sizeof(buckets));

        /* Assign primitives to buckets */
        float centroid_min = (split_axis == 0) ? centroid_bounds.min.x :
                            (split_axis == 1) ? centroid_bounds.min.y : centroid_bounds.min.z;
        float centroid_max = (split_axis == 0) ? centroid_bounds.max.x :
                            (split_axis == 1) ? centroid_bounds.max.y : centroid_bounds.max.z;
        float extent = centroid_max - centroid_min;

        if (extent > 1e-6f) {
            for (uint32_t i = start; i < end; i++) {
                float centroid_val = (split_axis == 0) ? primitives[i].centroid.x :
                                    (split_axis == 1) ? primitives[i].centroid.y :
                                                        primitives[i].centroid.z;
                uint32_t bucket_idx = (uint32_t)(GEOMLIB_BVH_SAH_BUCKETS *
                    ((centroid_val - centroid_min) / extent));
                if (bucket_idx >= GEOMLIB_BVH_SAH_BUCKETS) bucket_idx = GEOMLIB_BVH_SAH_BUCKETS - 1;

                buckets[bucket_idx].count++;
                if (buckets[bucket_idx].count == 1) {
                    buckets[bucket_idx].bounds = primitives[i].bounds;
                } else {
                    buckets[bucket_idx].bounds = aabb_union(&buckets[bucket_idx].bounds, &primitives[i].bounds);
                }
            }

            /* Compute costs for splitting at each bucket boundary */
            float min_cost = FLT_MAX;
            uint32_t min_cost_split = 0;

            for (uint32_t i = 0; i < GEOMLIB_BVH_SAH_BUCKETS - 1; i++) {
                geomlib_aabb left_bounds = buckets[0].bounds;
                uint32_t left_count = buckets[0].count;
                for (uint32_t j = 1; j <= i; j++) {
                    if (buckets[j].count > 0) {
                        left_bounds = aabb_union(&left_bounds, &buckets[j].bounds);
                        left_count += buckets[j].count;
                    }
                }

                geomlib_aabb right_bounds = buckets[i + 1].bounds;
                uint32_t right_count = buckets[i + 1].count;
                for (uint32_t j = i + 2; j < GEOMLIB_BVH_SAH_BUCKETS; j++) {
                    if (buckets[j].count > 0) {
                        right_bounds = aabb_union(&right_bounds, &buckets[j].bounds);
                        right_count += buckets[j].count;
                    }
                }

                if (left_count > 0 && right_count > 0) {
                    float cost = GEOMLIB_BVH_TRAVERSAL_COST +
                                (left_count * aabb_surface_area(&left_bounds) +
                                 right_count * aabb_surface_area(&right_bounds)) /
                                aabb_surface_area(&bounds) * GEOMLIB_BVH_INTERSECTION_COST;

                    if (cost < min_cost) {
                        min_cost = cost;
                        min_cost_split = i;
                    }
                }
            }

            /* Partition primitives based on best split */
            bvh_primitive_info *pmid = primitives + start;
            for (uint32_t i = start; i < end; i++) {
                float centroid_val = (split_axis == 0) ? primitives[i].centroid.x :
                                    (split_axis == 1) ? primitives[i].centroid.y :
                                                        primitives[i].centroid.z;
                uint32_t bucket_idx = (uint32_t)(GEOMLIB_BVH_SAH_BUCKETS *
                    ((centroid_val - centroid_min) / extent));
                if (bucket_idx >= GEOMLIB_BVH_SAH_BUCKETS) bucket_idx = GEOMLIB_BVH_SAH_BUCKETS - 1;

                if (bucket_idx <= min_cost_split) {
                    bvh_primitive_info tmp = *pmid;
                    *pmid = primitives[i];
                    primitives[i] = tmp;
                    pmid++;
                }
            }

            mid = (uint32_t)(pmid - primitives);
            if (mid == start || mid == end) mid = start + prim_count / 2;
        }
    } else {
        /* Median split (fast) */
        if (split_axis == 0) {
            qsort(primitives + start, prim_count, sizeof(bvh_primitive_info), compare_primitive_x);
        } else if (split_axis == 1) {
            qsort(primitives + start, prim_count, sizeof(bvh_primitive_info), compare_primitive_y);
        } else {
            qsort(primitives + start, prim_count, sizeof(bvh_primitive_info), compare_primitive_z);
        }
    }

    /* Recursively build children */
    uint32_t left = build_bvh_recursive(bvh, primitives, start, mid, quality);
    uint32_t right = build_bvh_recursive(bvh, primitives, mid, end, quality);

    node = &bvh->nodes[node_idx];  /* Refresh pointer after recursion */
    node->data.interior.left_child = left;
    node->data.interior.right_child = right;

    return node_idx;
}

/* ============================================================================
 * BVH Construction
 * ============================================================================ */

int geomlib_bvh_build(
    geomlib_bvh *bvh,
    const geomlib_trimesh *mesh,
    geomlib_bvh_build_quality quality
) {
    if (!bvh || !mesh || mesh->triangle_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    bvh->mesh = mesh;
    bvh->primitive_count = mesh->triangle_count;
    bvh->node_count = 0;
    bvh->node_capacity = mesh->triangle_count * 2;  /* Estimate */

    bvh->nodes = (geomlib_bvh_node *)malloc(sizeof(geomlib_bvh_node) * bvh->node_capacity);
    bvh->primitive_indices = (uint32_t *)malloc(sizeof(uint32_t) * bvh->primitive_count);

    if (!bvh->nodes || !bvh->primitive_indices) {
        geomlib_bvh_free(bvh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Build primitive info array */
    bvh_primitive_info *primitives = (bvh_primitive_info *)malloc(
        sizeof(bvh_primitive_info) * mesh->triangle_count
    );

    if (!primitives) {
        geomlib_bvh_free(bvh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        primitives[i].primitive_index = (uint32_t)i;
        primitives[i].bounds = compute_triangle_bounds(mesh, (uint32_t)i);
        primitives[i].centroid = aabb_centroid(&primitives[i].bounds);
    }

    /* Build BVH */
    build_bvh_recursive(bvh, primitives, 0, (uint32_t)mesh->triangle_count, quality);

    /* Copy reordered primitive indices */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        bvh->primitive_indices[i] = primitives[i].primitive_index;
    }

    free(primitives);

    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_bvh_free(geomlib_bvh *bvh) {
    if (bvh) {
        free(bvh->nodes);
        free(bvh->primitive_indices);
        bvh->nodes = NULL;
        bvh->primitive_indices = NULL;
        bvh->node_count = 0;
        bvh->node_capacity = 0;
        bvh->primitive_count = 0;
        bvh->mesh = NULL;
    }
}

static void compute_bvh_stats_recursive(
    const geomlib_bvh *bvh,
    uint32_t node_idx,
    uint32_t depth,
    uint32_t *max_depth,
    uint32_t *leaf_count,
    uint64_t *total_leaf_prims
) {
    const geomlib_bvh_node *node = &bvh->nodes[node_idx];

    if (depth > *max_depth) {
        *max_depth = depth;
    }

    if (node->is_leaf) {
        (*leaf_count)++;
        *total_leaf_prims += node->data.leaf.primitive_count;
    } else {
        compute_bvh_stats_recursive(bvh, node->data.interior.left_child, depth + 1,
                                    max_depth, leaf_count, total_leaf_prims);
        compute_bvh_stats_recursive(bvh, node->data.interior.right_child, depth + 1,
                                    max_depth, leaf_count, total_leaf_prims);
    }
}

void geomlib_bvh_get_stats(
    const geomlib_bvh *bvh,
    uint32_t *max_depth_out,
    uint32_t *leaf_count_out,
    float *avg_leaf_size_out
) {
    uint32_t max_depth = 0;
    uint32_t leaf_count = 0;
    uint64_t total_leaf_prims = 0;

    if (bvh && bvh->node_count > 0) {
        compute_bvh_stats_recursive(bvh, 0, 0, &max_depth, &leaf_count, &total_leaf_prims);
    }

    if (max_depth_out) *max_depth_out = max_depth;
    if (leaf_count_out) *leaf_count_out = leaf_count;
    if (avg_leaf_size_out) {
        *avg_leaf_size_out = (leaf_count > 0) ? ((float)total_leaf_prims / leaf_count) : 0.0f;
    }
}

/* ============================================================================
 * BVH Ray Queries
 * ============================================================================ */

static bool traverse_bvh_ray(
    const geomlib_bvh *bvh,
    const geomlib_ray *ray,
    geomlib_ray_hit *closest_hit,
    bool any_hit
) {
    if (bvh->node_count == 0) {
        return false;
    }

    bool found_hit = false;
    uint32_t stack[64];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    geomlib_ray current_ray = *ray;

    while (stack_ptr > 0) {
        uint32_t node_idx = stack[--stack_ptr];
        const geomlib_bvh_node *node = &bvh->nodes[node_idx];

        if (!geomlib_ray_aabb_intersect(&current_ray, &node->bounds)) {
            continue;
        }

        if (node->is_leaf) {
            /* Test primitives in leaf */
            for (uint32_t i = 0; i < node->data.leaf.primitive_count; i++) {
                uint32_t prim_idx = bvh->primitive_indices[node->data.leaf.primitive_offset + i];

                uint32_t i0 = bvh->mesh->indices[prim_idx * 3 + 0];
                uint32_t i1 = bvh->mesh->indices[prim_idx * 3 + 1];
                uint32_t i2 = bvh->mesh->indices[prim_idx * 3 + 2];

                geomlib_vec3 v0 = bvh->mesh->vertices[i0];
                geomlib_vec3 v1 = bvh->mesh->vertices[i1];
                geomlib_vec3 v2 = bvh->mesh->vertices[i2];

                float t, u, v;
                if (geomlib_ray_triangle_intersect(&current_ray, v0, v1, v2, &t, &u, &v)) {
                    if (any_hit) {
                        return true;
                    }

                    found_hit = true;
                    current_ray.t_max = t;

                    if (closest_hit) {
                        closest_hit->hit = true;
                        closest_hit->t = t;
                        closest_hit->position = geomlib_ray_at(&current_ray, t);
                        closest_hit->barycentric.x = 1.0f - u - v;
                        closest_hit->barycentric.y = u;
                        closest_hit->barycentric.z = v;
                        closest_hit->primitive_id = prim_idx;

                        /* Interpolate normal */
                        if (bvh->mesh->normals) {
                            geomlib_vec3 n0 = bvh->mesh->normals[i0];
                            geomlib_vec3 n1 = bvh->mesh->normals[i1];
                            geomlib_vec3 n2 = bvh->mesh->normals[i2];

                            closest_hit->normal.x = closest_hit->barycentric.x * n0.x +
                                                   closest_hit->barycentric.y * n1.x +
                                                   closest_hit->barycentric.z * n2.x;
                            closest_hit->normal.y = closest_hit->barycentric.x * n0.y +
                                                   closest_hit->barycentric.y * n1.y +
                                                   closest_hit->barycentric.z * n2.y;
                            closest_hit->normal.z = closest_hit->barycentric.x * n0.z +
                                                   closest_hit->barycentric.y * n1.z +
                                                   closest_hit->barycentric.z * n2.z;

                            closest_hit->normal = geomlib_vec3_normalize(closest_hit->normal);
                        } else {
                            /* Compute face normal */
                            geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
                            geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
                            closest_hit->normal = geomlib_vec3_normalize(geomlib_vec3_cross(e1, e2));
                        }
                    }
                }
            }
        } else {
            /* Push children to stack */
            if (stack_ptr < 63) {
                stack[stack_ptr++] = node->data.interior.left_child;
                stack[stack_ptr++] = node->data.interior.right_child;
            }
        }
    }

    return found_hit;
}

bool geomlib_bvh_trace_ray(
    const geomlib_bvh *bvh,
    const geomlib_ray *ray,
    geomlib_ray_hit *hit
) {
    if (hit) {
        hit->hit = false;
    }
    return traverse_bvh_ray(bvh, ray, hit, false);
}

bool geomlib_bvh_test_occlusion(
    const geomlib_bvh *bvh,
    const geomlib_ray *ray
) {
    return traverse_bvh_ray(bvh, ray, NULL, true);
}

/* ============================================================================
 * BVH Point Queries
 * ============================================================================ */

bool geomlib_bvh_closest_point(
    const geomlib_bvh *bvh,
    geomlib_vec3 query,
    geomlib_vec3 *closest_out,
    float *distance_out,
    uint32_t *primitive_id_out
) {
    if (!bvh || bvh->node_count == 0 || !closest_out) {
        return false;
    }

    geomlib_vec3 closest_point = {0, 0, 0};
    float min_dist_sq = FLT_MAX;
    uint32_t closest_prim = 0;

    uint32_t stack[64];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        uint32_t node_idx = stack[--stack_ptr];
        const geomlib_bvh_node *node = &bvh->nodes[node_idx];

        /* Check if node could contain closer point */
        geomlib_vec3 node_closest = geomlib_aabb_closest_point(&node->bounds, query);
        geomlib_vec3 diff = geomlib_vec3_sub(node_closest, query);
        float node_dist_sq = geomlib_vec3_dot(diff, diff);

        if (node_dist_sq >= min_dist_sq) {
            continue;
        }

        if (node->is_leaf) {
            /* Test primitives in leaf */
            for (uint32_t i = 0; i < node->data.leaf.primitive_count; i++) {
                uint32_t prim_idx = bvh->primitive_indices[node->data.leaf.primitive_offset + i];

                uint32_t i0 = bvh->mesh->indices[prim_idx * 3 + 0];
                uint32_t i1 = bvh->mesh->indices[prim_idx * 3 + 1];
                uint32_t i2 = bvh->mesh->indices[prim_idx * 3 + 2];

                geomlib_vec3 v0 = bvh->mesh->vertices[i0];
                geomlib_vec3 v1 = bvh->mesh->vertices[i1];
                geomlib_vec3 v2 = bvh->mesh->vertices[i2];

                geomlib_vec3 pt = geomlib_closest_point_triangle(query, v0, v1, v2);
                geomlib_vec3 diff_pt = geomlib_vec3_sub(pt, query);
                float dist_sq = geomlib_vec3_dot(diff_pt, diff_pt);

                if (dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                    closest_point = pt;
                    closest_prim = prim_idx;
                }
            }
        } else {
            /* Push children to stack (closer one last for DFS) */
            const geomlib_bvh_node *left = &bvh->nodes[node->data.interior.left_child];
            const geomlib_bvh_node *right = &bvh->nodes[node->data.interior.right_child];

            geomlib_vec3 left_closest = geomlib_aabb_closest_point(&left->bounds, query);
            geomlib_vec3 right_closest = geomlib_aabb_closest_point(&right->bounds, query);

            geomlib_vec3 left_diff = geomlib_vec3_sub(left_closest, query);
            geomlib_vec3 right_diff = geomlib_vec3_sub(right_closest, query);

            float left_dist_sq = geomlib_vec3_dot(left_diff, left_diff);
            float right_dist_sq = geomlib_vec3_dot(right_diff, right_diff);

            if (left_dist_sq < right_dist_sq) {
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.right_child;
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.left_child;
            } else {
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.left_child;
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.right_child;
            }
        }
    }

    *closest_out = closest_point;
    if (distance_out) *distance_out = sqrtf(min_dist_sq);
    if (primitive_id_out) *primitive_id_out = closest_prim;

    return true;
}

typedef struct knn_candidate {
    float dist2;
    geomlib_vec3 point;
    uint32_t primitive;
} knn_candidate;

static int compare_knn_candidates(const void *left, const void *right) {
    const knn_candidate *a = left, *b = right;
    return (a->dist2 > b->dist2) - (a->dist2 < b->dist2);
}

size_t geomlib_bvh_knn_query(
    const geomlib_bvh *bvh,
    geomlib_vec3 query,
    size_t k,
    geomlib_vec3 *points_out,
    float *distances_out,
    uint32_t *primitive_ids_out
) {
    if (!bvh || !bvh->mesh || bvh->node_count == 0 || k == 0 || !points_out || k > SIZE_MAX / sizeof(knn_candidate))
        return 0;
    knn_candidate *best = malloc(k * sizeof(*best));
    if (!best) return 0;
    size_t count = 0;

    uint32_t stack[64];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;
    while (stack_ptr > 0) {
        uint32_t node_index = stack[--stack_ptr];
        const geomlib_bvh_node *node = &bvh->nodes[node_index];
        geomlib_vec3 node_point = geomlib_aabb_closest_point(&node->bounds, query);
        geomlib_vec3 node_delta = geomlib_vec3_sub(node_point, query);
        float node_dist2 = geomlib_vec3_dot(node_delta, node_delta);
        if (count == k) {
            float worst = 0.0f;
            for (size_t i = 0; i < count; i++) if (best[i].dist2 > worst) worst = best[i].dist2;
            if (node_dist2 > worst) continue;
        }

        if (node->is_leaf) {
            for (uint32_t i = 0; i < node->data.leaf.primitive_count; i++) {
                uint32_t primitive = bvh->primitive_indices[node->data.leaf.primitive_offset + i];
                uint32_t i0 = bvh->mesh->indices[primitive * 3 + 0];
                uint32_t i1 = bvh->mesh->indices[primitive * 3 + 1];
                uint32_t i2 = bvh->mesh->indices[primitive * 3 + 2];
                geomlib_vec3 point = geomlib_closest_point_triangle(query,
                    bvh->mesh->vertices[i0], bvh->mesh->vertices[i1],
                    bvh->mesh->vertices[i2]);
                geomlib_vec3 delta = geomlib_vec3_sub(point, query);
                float dist2 = geomlib_vec3_dot(delta, delta);
                bool duplicate = false;
                for (size_t j = 0; j < count && !duplicate; j++)
                    duplicate = best[j].primitive == primitive;
                if (duplicate) continue;
                if (count < k) {
                    best[count++] = (knn_candidate){dist2, point, primitive};
                    continue;
                }
                size_t worst = 0;
                for (size_t j = 1; j < count; j++)
                    if (best[j].dist2 > best[worst].dist2) worst = j;
                if (dist2 < best[worst].dist2)
                    best[worst] = (knn_candidate){dist2, point, primitive};
            }
        } else {
            const geomlib_bvh_node *left = &bvh->nodes[node->data.interior.left_child];
            const geomlib_bvh_node *right = &bvh->nodes[node->data.interior.right_child];
            geomlib_vec3 left_delta = geomlib_vec3_sub(
                geomlib_aabb_closest_point(&left->bounds, query), query);
            geomlib_vec3 right_delta = geomlib_vec3_sub(
                geomlib_aabb_closest_point(&right->bounds, query), query);
            float left_dist2 = geomlib_vec3_dot(left_delta, left_delta);
            float right_dist2 = geomlib_vec3_dot(right_delta, right_delta);
            /* Push the farther child first so the nearer subtree is visited next. */
            if (left_dist2 < right_dist2) {
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.right_child;
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.left_child;
            } else {
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.left_child;
                if (stack_ptr < 63) stack[stack_ptr++] = node->data.interior.right_child;
            }
        }
    }

    qsort(best, count, sizeof(*best), compare_knn_candidates);
    for (size_t i = 0; i < count; i++) {
        points_out[i] = best[i].point;
        if (distances_out) distances_out[i] = sqrtf(best[i].dist2);
        if (primitive_ids_out) primitive_ids_out[i] = best[i].primitive;
    }
    free(best);
    return count;
}

size_t geomlib_bvh_radius_query(
    const geomlib_bvh *bvh,
    geomlib_vec3 query,
    float radius,
    size_t max_results,
    geomlib_vec3 *points_out,
    float *distances_out,
    uint32_t *primitive_ids_out
) {
    if (!bvh || bvh->node_count == 0 || max_results == 0 || !points_out) {
        return 0;
    }

    float radius_sq = radius * radius;
    size_t found_count = 0;

    uint32_t stack[64];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0 && found_count < max_results) {
        uint32_t node_idx = stack[--stack_ptr];
        const geomlib_bvh_node *node = &bvh->nodes[node_idx];

        /* Check if node intersects query sphere */
        geomlib_vec3 node_closest = geomlib_aabb_closest_point(&node->bounds, query);
        geomlib_vec3 diff = geomlib_vec3_sub(node_closest, query);
        float node_dist_sq = geomlib_vec3_dot(diff, diff);

        if (node_dist_sq > radius_sq) {
            continue;
        }

        if (node->is_leaf) {
            /* Test primitives in leaf */
            for (uint32_t i = 0; i < node->data.leaf.primitive_count && found_count < max_results; i++) {
                uint32_t prim_idx = bvh->primitive_indices[node->data.leaf.primitive_offset + i];

                uint32_t i0 = bvh->mesh->indices[prim_idx * 3 + 0];
                uint32_t i1 = bvh->mesh->indices[prim_idx * 3 + 1];
                uint32_t i2 = bvh->mesh->indices[prim_idx * 3 + 2];

                geomlib_vec3 v0 = bvh->mesh->vertices[i0];
                geomlib_vec3 v1 = bvh->mesh->vertices[i1];
                geomlib_vec3 v2 = bvh->mesh->vertices[i2];

                geomlib_vec3 pt = geomlib_closest_point_triangle(query, v0, v1, v2);
                geomlib_vec3 diff_pt = geomlib_vec3_sub(pt, query);
                float dist_sq = geomlib_vec3_dot(diff_pt, diff_pt);

                if (dist_sq <= radius_sq) {
                    points_out[found_count] = pt;
                    if (distances_out) distances_out[found_count] = sqrtf(dist_sq);
                    if (primitive_ids_out) primitive_ids_out[found_count] = prim_idx;
                    found_count++;
                }
            }
        } else {
            /* Push children to stack */
            if (stack_ptr < 63) {
                stack[stack_ptr++] = node->data.interior.left_child;
                stack[stack_ptr++] = node->data.interior.right_child;
            }
        }
    }

    return found_count;
}
