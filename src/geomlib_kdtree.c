/*
 * geomlib - K-d tree (see geomlib_kdtree.h)
 *
 * Clean-room implementation: median-split recursive build on the widest
 * coordinate axis, with branch-and-bound kNN and radius queries.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint32_t index;
    int axis;
    int left, right;   /* child node indices, -1 for none */
    geomlib_vec3 split;
} kd_node;

struct geomlib_kdtree {
    const geomlib_vec3 *points;
    size_t count;
    kd_node *nodes;
    size_t node_count;
};

typedef struct { size_t index; size_t axis; } item;

static void sort_items(item *items, size_t n, int axis,
                       const geomlib_vec3 *points)
{
    size_t i, j;
    for (i = 1; i < n; ++i) {
        item cur = items[i];
        float key = axis==0 ? points[cur.index].x
                  : (axis==1 ? points[cur.index].y : points[cur.index].z);
        j = i;
        while (j > 0) {
            float pk = axis==0 ? points[items[j-1].index].x
                     : (axis==1 ? points[items[j-1].index].y
                                : points[items[j-1].index].z);
            if (pk <= key) break;
            items[j] = items[j-1];
            --j;
        }
        items[j] = cur;
    }
}

static int build_node(geomlib_kdtree *t, item *items, size_t n, int depth)
{
    int id;
    size_t mid;
    if (!n) return -1;
    if (t->node_count >= t->count*2+1) return -1;
    id = (int)t->node_count++;
    {
        geomlib_vec3 lo = t->points[items[0].index];
        geomlib_vec3 hi=lo;
        size_t i;
        for (i = 1; i < n; ++i) {
            geomlib_vec3 p = t->points[items[i].index];
            if (p.x<lo.x) lo.x=p.x;
            if (p.y<lo.y) lo.y=p.y;
            if (p.z<lo.z) lo.z=p.z;
            if (p.x>hi.x) hi.x=p.x;
            if (p.y>hi.y) hi.y=p.y;
            if (p.z>hi.z) hi.z=p.z;
        }
        {
            geomlib_vec3 s = geomlib_vec3_create(hi.x-lo.x, hi.y-lo.y, hi.z-lo.z);
            int axis = s.x > s.y ? (s.x > s.z ? 0 : 2) : (s.y > s.z ? 1 : 2);
            sort_items(items, n, axis, t->points);
            mid = n/2;
            t->nodes[id].axis = axis;
            t->nodes[id].index = items[mid].index;
            t->nodes[id].split = t->points[items[mid].index];
        }
    }
    t->nodes[id].left = build_node(t, items, mid, depth+1);
    t->nodes[id].right = build_node(t, items+mid+1, n-mid-1, depth+1);
    return id;
}

geomlib_kdtree *geomlib_kdtree_build(const geomlib_vec3 *points, size_t count)
{
    geomlib_kdtree *t;
    item *items;
    size_t i;
    if (!points || !count) return NULL;
    t = (geomlib_kdtree *)calloc(1, sizeof(*t));
    if (!t) return NULL;
    t->points = points;
    t->count = count;
    t->nodes = (kd_node *)calloc(count*2+1, sizeof(kd_node));
    items = (item *)malloc(count*sizeof(item));
    if (!t->nodes || !items) { free(items); geomlib_kdtree_free(t); return NULL; }
    for (i = 0; i < count; ++i) { items[i].index = i; items[i].axis = 0; }
    build_node(t, items, count, 0);
    free(items);
    return t;
}

void geomlib_kdtree_free(geomlib_kdtree *t)
{
    if (!t) return;
    free(t->nodes);
    free(t);
}

static float sqdist(geomlib_vec3 a, geomlib_vec3 b)
{
    geomlib_vec3 d = geomlib_vec3_sub(a,b);
    return geomlib_vec3_dot(d,d);
}

static void knn_walk(const geomlib_kdtree *t, int id, geomlib_vec3 q,
                     size_t k, size_t *found, uint32_t *idx, float *dst)
{
    const kd_node *n;
    float d, worst, plane;
    int near, far;
    size_t i;
    if (id < 0) return;
    n = &t->nodes[id];
    d = sqdist(q, n->split);
    if (*found < k || d < dst[k-1]) {
        size_t pos = *found < k ? (*found)++ : k-1;
        idx[pos] = n->index;
        dst[pos] = d;
        for (i = pos; i > 0 && dst[i] < dst[i-1]; --i) {
            float td=dst[i]; uint32_t ti=idx[i];
            dst[i]=dst[i-1]; idx[i]=idx[i-1];
            dst[i-1]=td; idx[i-1]=ti;
        }
    }
    worst = *found < k ? INFINITY : dst[k-1];
    plane = n->axis==0 ? q.x-n->split.x : (n->axis==1 ? q.y-n->split.y : q.z-n->split.z);
    near = plane <= 0 ? n->left : n->right;
    far = plane <= 0 ? n->right : n->left;
    knn_walk(t, near, q, k, found, idx, dst);
    worst = *found < k ? INFINITY : dst[k-1];
    if (plane*plane <= worst) knn_walk(t, far, q, k, found, idx, dst);
}

size_t geomlib_kdtree_knn(const geomlib_kdtree *t, geomlib_vec3 q,
                          size_t k, uint32_t *indices, float *distances)
{
    size_t found = 0, i;
    float *tmp;
    if (!t || !k || !indices) return 0;
    tmp = distances ? distances : (float *)malloc(k*sizeof(float));
    if (!tmp) return 0;
    for (i = 0; i < k; ++i) tmp[i] = INFINITY;
    knn_walk(t, 0, q, k, &found, indices, tmp);
    if (distances) {
        for (i = 0; i < found; ++i) distances[i] = sqrtf(tmp[i]);
    } else {
        free(tmp);
    }
    return found;
}

static void radius_walk(const geomlib_kdtree *t, int id, geomlib_vec3 q,
                        float r2, uint32_t *idx, float *dst,
                        size_t max_results, size_t *found)
{
    const kd_node *n;
    float d, plane;
    if (id < 0 || *found >= max_results) return;
    n = &t->nodes[id];
    d = sqdist(q, n->split);
    if (d <= r2 && *found < max_results) {
        idx[*found] = n->index;
        if (dst) dst[*found] = sqrtf(d);
        (*found)++;
    }
    plane = n->axis==0 ? q.x-n->split.x : (n->axis==1 ? q.y-n->split.y : q.z-n->split.z);
    radius_walk(t, plane <= 0 ? n->left : n->right, q, r2, idx, dst, max_results, found);
    if (plane*plane <= r2)
        radius_walk(t, plane <= 0 ? n->right : n->left, q, r2, idx, dst, max_results, found);
}

size_t geomlib_kdtree_radius(const geomlib_kdtree *t, geomlib_vec3 q,
                             float radius, uint32_t *indices,
                             float *distances, size_t max_results)
{
    size_t found = 0;
    if (!t || !indices || !(radius > 0.0f)) return 0;
    radius_walk(t, 0, q, radius*radius, indices, distances, max_results, &found);
    return found;
}
