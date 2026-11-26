/*
 * geomlib - Point-cloud octree (see geomlib_octree.h)
 *
 * Clean-room implementation: recursively split nodes with more than
 * `leaf_capacity` points into eight octants, then prune radius/kNN queries by
 * the distance from the query point to each node's bounding box.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_octree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    geomlib_aabb bounds;
    int children[8];
    uint32_t *points;
    size_t count;
    int depth;
} oct_node;

struct geomlib_octree {
    const geomlib_vec3 *points;
    size_t count;
    oct_node *nodes;
    size_t node_count, node_capacity;
    geomlib_octree_params params;
};

geomlib_octree_params geomlib_octree_params_default(void)
{
    geomlib_octree_params p;
    p.leaf_capacity = 8;
    p.max_depth = 12;
    return p;
}

static int child_index(geomlib_aabb box, geomlib_vec3 p)
{
    geomlib_vec3 c = geomlib_aabb_center(box);
    return (p.x >= c.x ? 1 : 0) | (p.y >= c.y ? 2 : 0) | (p.z >= c.z ? 4 : 0);
}

static geomlib_aabb child_bounds(geomlib_aabb box, int oct)
{
    geomlib_vec3 c = geomlib_aabb_center(box);
    geomlib_vec3 lo = box.min, hi = box.max;
    if (oct & 1) lo.x = c.x; else hi.x = c.x;
    if (oct & 2) lo.y = c.y; else hi.y = c.y;
    if (oct & 4) lo.z = c.z; else hi.z = c.z;
    return geomlib_aabb_create(lo, hi);
}

static int add_node(geomlib_octree *t)
{
    if (t->node_count == t->node_capacity) {
        size_t nc = t->node_capacity ? t->node_capacity*2 : 64;
        oct_node *nn = (oct_node *)realloc(t->nodes, nc*sizeof(oct_node));
        if (!nn) return -1;
        t->nodes = nn;
        t->node_capacity = nc;
    }
    memset(&t->nodes[t->node_count], 0, sizeof(oct_node));
    for (int i = 0; i < 8; ++i) t->nodes[t->node_count].children[i] = -1;
    return (int)t->node_count++;
}

static void split_node(geomlib_octree *t, int id)
{
    oct_node *n = &t->nodes[id];
    size_t counts[8] = {0};
    uint32_t *assign = (uint32_t *)malloc(n->count*sizeof(uint32_t));
    if (!assign) return;
    for (size_t i = 0; i < n->count; ++i) {
        assign[i] = (uint32_t)child_index(n->bounds, t->points[n->points[i]]);
        counts[assign[i]]++;
    }
    for (int c = 0; c < 8; ++c) {
        size_t pos[8] = {0}, j;
        int child;
        if (!counts[c]) continue;
        child = add_node(t);
        if (child < 0) break;
        n = &t->nodes[id];
        n->children[c] = child;
        t->nodes[child].bounds = child_bounds(n->bounds, c);
        t->nodes[child].depth = n->depth+1;
        t->nodes[child].points = (uint32_t *)malloc(counts[c]*sizeof(uint32_t));
        if (!t->nodes[child].points) { t->nodes[child].count = 0; continue; }
        t->nodes[child].count = counts[c];
        for (j = 0; j < n->count; ++j) {
            if (assign[j] != (uint32_t)c) continue;
            t->nodes[child].points[pos[c]++] = n->points[j];
        }
    }
    free(assign);
    n = &t->nodes[id];
    free(n->points);
    n->points = NULL;
    n->count = 0;
    /* Cache child ids first: recursive splits may realloc the node array. */
    {
        int children[8];
        for (int c = 0; c < 8; ++c) children[c] = n->children[c];
        for (int c = 0; c < 8; ++c) {
            if (children[c] < 0) continue;
            if (t->nodes[children[c]].count > t->params.leaf_capacity &&
                t->nodes[children[c]].depth < t->params.max_depth)
                split_node(t, children[c]);
        }
    }
}

geomlib_octree *geomlib_octree_build(const geomlib_vec3 *points, size_t count,
                                     const geomlib_octree_params *params)
{
    geomlib_octree *t;
    int root;
    if (!points || !count) return NULL;
    t = (geomlib_octree *)calloc(1, sizeof(*t));
    if (!t) return NULL;
    t->points = points;
    t->count = count;
    t->params = params ? *params : geomlib_octree_params_default();
    if (t->params.leaf_capacity < 1) t->params.leaf_capacity = 1;
    if (t->params.max_depth < 1) t->params.max_depth = 1;
    if (t->params.max_depth > 20) t->params.max_depth = 20;
    root = add_node(t);
    if (root < 0) { geomlib_octree_free(t); return NULL; }
    {
        geomlib_vec3 lo=points[0], hi=points[0];
        size_t i;
        for (i = 1; i < count; ++i) {
            if (points[i].x<lo.x) lo.x=points[i].x;
            if (points[i].y<lo.y) lo.y=points[i].y;
            if (points[i].z<lo.z) lo.z=points[i].z;
            if (points[i].x>hi.x) hi.x=points[i].x;
            if (points[i].y>hi.y) hi.y=points[i].y;
            if (points[i].z>hi.z) hi.z=points[i].z;
        }
        {
            geomlib_vec3 pad = geomlib_vec3_mul(
                geomlib_vec3_sub(hi,lo), 1e-3f);
            if (geomlib_vec3_length(pad) < 1e-9f) pad = geomlib_vec3_create(1e-3f,1e-3f,1e-3f);
            t->nodes[root].bounds = geomlib_aabb_create(
                geomlib_vec3_sub(lo,pad), geomlib_vec3_add(hi,pad));
        }
    }
    t->nodes[root].points = (uint32_t *)malloc(count*sizeof(uint32_t));
    if (!t->nodes[root].points) { geomlib_octree_free(t); return NULL; }
    t->nodes[root].count = count;
    for (size_t i = 0; i < count; ++i) t->nodes[root].points[i] = (uint32_t)i;
    if (count > t->params.leaf_capacity) split_node(t, root);
    return t;
}

void geomlib_octree_free(geomlib_octree *t)
{
    size_t i;
    if (!t) return;
    for (i = 0; i < t->node_count; ++i) free(t->nodes[i].points);
    free(t->nodes);
    free(t);
}

static float box_dist2(geomlib_aabb b, geomlib_vec3 q)
{
    float dx = fmaxf(fmaxf(b.min.x-q.x, q.x-b.max.x), 0.0f);
    float dy = fmaxf(fmaxf(b.min.y-q.y, q.y-b.max.y), 0.0f);
    float dz = fmaxf(fmaxf(b.min.z-q.z, q.z-b.max.z), 0.0f);
    return dx*dx+dy*dy+dz*dz;
}

static void radius_walk(const geomlib_octree *t, int id, geomlib_vec3 q,
                        float r2, uint32_t *idx, float *dst,
                        size_t max_results, size_t *found)
{
    const oct_node *n;
    size_t i;
    if (id < 0 || *found >= max_results) return;
    n = &t->nodes[id];
    if (box_dist2(n->bounds, q) > r2) return;
    for (i = 0; i < n->count && *found < max_results; ++i) {
        uint32_t p = n->points[i];
        geomlib_vec3 d = geomlib_vec3_sub(t->points[p], q);
        float dd = geomlib_vec3_dot(d,d);
        if (dd <= r2) {
            idx[*found] = p;
            if (dst) dst[*found] = sqrtf(dd);
            (*found)++;
        }
    }
    for (int c = 0; c < 8 && *found < max_results; ++c)
        radius_walk(t, n->children[c], q, r2, idx, dst, max_results, found);
}

size_t geomlib_octree_radius(const geomlib_octree *t, geomlib_vec3 q,
                             float radius, uint32_t *indices,
                             float *distances, size_t max_results)
{
    size_t found = 0;
    if (!t || !indices || !(radius > 0.0f)) return 0;
    radius_walk(t, 0, q, radius*radius, indices, distances, max_results, &found);
    return found;
}

typedef struct { int id; float d2; } node_dist;

static int cmp_node_dist(const void *a, const void *b)
{
    const node_dist *x=(const node_dist*)a, *y=(const node_dist*)b;
    return x->d2 < y->d2 ? -1 : (x->d2 > y->d2 ? 1 : 0);
}

static void knn_walk(const geomlib_octree *t, geomlib_vec3 q, size_t k,
                     uint32_t *idx, float *dst, size_t *found)
{
    node_dist *queue = (node_dist *)malloc((t->node_count+1)*sizeof(node_dist));
    size_t head=0, tail=0;
    if (!queue) return;
    queue[tail].id = 0;
    queue[tail++].d2 = box_dist2(t->nodes[0].bounds, q);
    while (head < tail) {
        const oct_node *n;
        size_t i;
        qsort(queue+head, tail-head, sizeof(node_dist), cmp_node_dist);
        n = &t->nodes[queue[head].id];
        ++head;
        for (i = 0; i < n->count; ++i) {
            uint32_t p = n->points[i];
            geomlib_vec3 d = geomlib_vec3_sub(t->points[p], q);
            float dd = geomlib_vec3_dot(d,d);
            size_t pos;
            if (*found < k || dd < dst[k-1]) {
                pos = *found < k ? (*found)++ : k-1;
                idx[pos] = p;
                dst[pos] = dd;
                while (pos > 0 && dst[pos] < dst[pos-1]) {
                    float td=dst[pos]; uint32_t ti=idx[pos];
                    dst[pos]=dst[pos-1]; idx[pos]=idx[pos-1];
                    dst[pos-1]=td; idx[pos-1]=ti;
                    --pos;
                }
            }
        }
        if (*found == k && queue[head < tail ? head : tail-1].d2 > dst[k-1]) break;
        for (int c = 0; c < 8; ++c) {
            int child = n->children[c];
            if (child < 0) continue;
            queue[tail].id = child;
            queue[tail++].d2 = box_dist2(t->nodes[child].bounds, q);
        }
    }
    free(queue);
}

size_t geomlib_octree_knn(const geomlib_octree *t, geomlib_vec3 q,
                          size_t k, uint32_t *indices, float *distances)
{
    size_t found = 0, i;
    float *tmp;
    if (!t || !k || !indices) return 0;
    tmp = distances ? distances : (float *)malloc(k*sizeof(float));
    if (!tmp) return 0;
    for (i = 0; i < k; ++i) tmp[i] = INFINITY;
    knn_walk(t, q, k, indices, tmp, &found);
    if (distances) for (i = 0; i < found; ++i) distances[i] = sqrtf(tmp[i]);
    else free(tmp);
    return found;
}
