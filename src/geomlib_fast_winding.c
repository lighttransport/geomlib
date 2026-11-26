/*
 * geomlib - Fast winding number (see geomlib_fast_winding.h)
 *
 * Clean-room implementation of the hierarchical far-field idea from Barill et
 * al., "Fast Winding Numbers for Soups and Clouds", SIGGRAPH 2018: cluster
 * triangles into a tree; far clusters are approximated by their normal
 * moment's dipole solid angle, near leaves are evaluated exactly.
 *
 * The tree holds a non-owning pointer to the source mesh, which must outlive
 * the tree. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_fast_winding.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI_F 3.14159265358979323846f

typedef struct {
    uint32_t begin, end;
    int left, right;
    geomlib_vec3 centroid;
    geomlib_vec3 moment;
    float radius;
} wn_node;

struct geomlib_winding_tree {
    const geomlib_trimesh *mesh;
    uint32_t *indices;
    wn_node *nodes;
    size_t node_count;
};

geomlib_winding_tree *geomlib_winding_tree_build(const geomlib_trimesh *mesh,
                                                 size_t leaf_size)
{
    geomlib_winding_tree *t;
    size_t max_nodes, count = 1, i;
    int *stack;
    size_t sp = 0;

    if (!mesh || !mesh->vertices || !mesh->indices || mesh->triangle_count == 0)
        return NULL;
    if (leaf_size == 0) leaf_size = 16;
    if (leaf_size > 512) leaf_size = 512;
    max_nodes = 4*mesh->triangle_count/leaf_size + 64;

    t = (geomlib_winding_tree *)calloc(1, sizeof(*t));
    if (!t) return NULL;
    t->mesh = mesh;
    t->indices = (uint32_t *)malloc(mesh->triangle_count*sizeof(uint32_t));
    t->nodes = (wn_node *)calloc(max_nodes, sizeof(wn_node));
    stack = (int *)malloc(max_nodes*sizeof(int));
    if (!t->indices || !t->nodes || !stack) {
        free(stack); free(t->nodes); free(t->indices); free(t);
        return NULL;
    }
    for (i = 0; i < mesh->triangle_count; ++i) t->indices[i] = (uint32_t)i;
    t->nodes[0].begin = 0;
    t->nodes[0].end = (uint32_t)mesh->triangle_count;
    t->nodes[0].left = t->nodes[0].right = -1;
    stack[sp++] = 0;

    while (sp) {
        int id = (int)stack[--sp];
        wn_node *n = &t->nodes[id];
        size_t len = n->end - n->begin;
        if (len <= leaf_size || count+2 >= max_nodes) continue;
        {
            geomlib_vec3 lo = {FLT_MAX,FLT_MAX,FLT_MAX}, hi = {-FLT_MAX,-FLT_MAX,-FLT_MAX};
            geomlib_vec3 size;
            int axis;
            size_t j, mid = n->begin + len/2;
            int cl = (int)count++, cr = (int)count++;
            for (j = n->begin; j < n->end; ++j) {
                geomlib_vec3 c = mesh->vertices[mesh->indices[t->indices[j]*3]];
                if (c.x<lo.x) lo.x=c.x;
                if (c.y<lo.y) lo.y=c.y;
                if (c.z<lo.z) lo.z=c.z;
                if (c.x>hi.x) hi.x=c.x;
                if (c.y>hi.y) hi.y=c.y;
                if (c.z>hi.z) hi.z=c.z;
            }
            size = geomlib_vec3_create(hi.x-lo.x, hi.y-lo.y, hi.z-lo.z);
            axis = size.x > size.y ? (size.x > size.z ? 0 : 2)
                                   : (size.y > size.z ? 1 : 2);
            /* Insertion sort by the chosen axis. */
            for (j = n->begin+1; j < n->end; ++j) {
                uint32_t cur = t->indices[j];
                geomlib_vec3 pc = mesh->vertices[mesh->indices[cur*3]];
                float key = axis == 0 ? pc.x : (axis == 1 ? pc.y : pc.z);
                size_t p = j;
                while (p > n->begin) {
                    uint32_t prev = t->indices[p-1];
                    geomlib_vec3 pp = mesh->vertices[mesh->indices[prev*3]];
                    float pk = axis == 0 ? pp.x : (axis == 1 ? pp.y : pp.z);
                    if (pk <= key) break;
                    t->indices[p] = prev;
                    --p;
                }
                t->indices[p] = cur;
            }
            t->nodes[cl].begin = n->begin;
            t->nodes[cl].end = (uint32_t)mid;
            t->nodes[cr].begin = (uint32_t)mid;
            t->nodes[cr].end = n->end;
            t->nodes[cl].left = t->nodes[cl].right = -1;
            t->nodes[cr].left = t->nodes[cr].right = -1;
            n->left = cl;
            n->right = cr;
            stack[sp++] = cl;
            stack[sp++] = cr;
        }
    }
    free(stack);
    t->node_count = count;

    for (i = 0; i < count; ++i) {
        wn_node *n = &t->nodes[i];
        geomlib_vec3 c = {0,0,0}, m = {0,0,0};
        float area = 0.0f;
        size_t j, k;
        for (j = n->begin; j < n->end; ++j) {
            uint32_t tri = t->indices[j];
            uint32_t a = mesh->indices[tri*3], b = mesh->indices[tri*3+1],
                     c2 = mesh->indices[tri*3+2];
            geomlib_vec3 e1 = geomlib_vec3_sub(mesh->vertices[b], mesh->vertices[a]);
            geomlib_vec3 e2 = geomlib_vec3_sub(mesh->vertices[c2], mesh->vertices[a]);
            geomlib_vec3 cr = geomlib_vec3_cross(e1, e2);
            float a2 = 0.5f*geomlib_vec3_length(cr);
            geomlib_vec3 tc = geomlib_vec3_mul(
                geomlib_vec3_add(geomlib_vec3_add(mesh->vertices[a], mesh->vertices[b]),
                                 mesh->vertices[c2]), 1.0f/3.0f);
            c = geomlib_vec3_add(c, geomlib_vec3_mul(tc, a2));
            m = geomlib_vec3_add(m, geomlib_vec3_mul(cr, 0.5f));
            area += a2;
        }
        n->centroid = area > 1e-12f ? geomlib_vec3_mul(c, 1.0f/area)
                                    : mesh->vertices[mesh->indices[t->indices[n->begin]*3]];
        n->moment = m;
        n->radius = 0.0f;
        for (j = n->begin; j < n->end; ++j) {
            uint32_t tri = t->indices[j];
            for (k = 0; k < 3; ++k) {
                geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[mesh->indices[tri*3+k]],
                                                  n->centroid);
                float dd = geomlib_vec3_length(d);
                if (dd > n->radius) n->radius = dd;
            }
        }
    }
    return t;
}

void geomlib_winding_tree_free(geomlib_winding_tree *t)
{
    if (!t) return;
    free(t->indices);
    free(t->nodes);
    free(t);
}

static float exact_winding(const geomlib_winding_tree *t, const wn_node *n,
                           geomlib_vec3 q)
{
    double sum = 0.0;
    size_t i;
    const geomlib_trimesh *mesh = t->mesh;
    for (i = n->begin; i < n->end; ++i) {
        uint32_t tri = t->indices[i];
        geomlib_vec3 a = geomlib_vec3_sub(mesh->vertices[mesh->indices[tri*3]], q);
        geomlib_vec3 b = geomlib_vec3_sub(mesh->vertices[mesh->indices[tri*3+1]], q);
        geomlib_vec3 c = geomlib_vec3_sub(mesh->vertices[mesh->indices[tri*3+2]], q);
        double la = sqrt((double)geomlib_vec3_dot(a,a));
        double lb = sqrt((double)geomlib_vec3_dot(b,b));
        double lc = sqrt((double)geomlib_vec3_dot(c,c));
        double num = geomlib_vec3_dot(a, geomlib_vec3_cross(b,c));
        double den = la*lb*lc + geomlib_vec3_dot(a,b)*lc +
                     geomlib_vec3_dot(b,c)*la + geomlib_vec3_dot(c,a)*lb;
        if (fabs(den) > 1e-18) sum += 2.0*atan2(num, den);
    }
    return (float)(sum/(4.0*PI_F));
}

static float walk(const geomlib_winding_tree *t, int id, geomlib_vec3 q)
{
    const wn_node *n = &t->nodes[id];
    geomlib_vec3 r = geomlib_vec3_sub(n->centroid, q);
    float dist = geomlib_vec3_length(r);
    double far;
    if (dist < 1e-8f) dist = 1e-8f;
    far = geomlib_vec3_dot(n->moment, r) / (4.0*PI_F*dist*dist*dist);
    if (n->left < 0) return exact_winding(t, n, q);
    if (dist > 5.0f*n->radius + 1e-6f) return (float)far;
    return walk(t, n->left, q) + walk(t, n->right, q);
}

float geomlib_fast_winding_number(const geomlib_winding_tree *t,
                                  geomlib_vec3 point)
{
    if (!t) return 0.0f;
    return walk(t, 0, point);
}
