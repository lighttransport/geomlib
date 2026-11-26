/*
 * geomlib - Dynamic BVH refit (see geomlib_bvh_refit.h)
 *
 * Clean-room implementation: recompute each leaf's AABB from its triangles in
 * the current mesh, then union child bounds bottom-up for interior nodes.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_bvh_refit.h"

#include <float.h>

static void refit_node(geomlib_bvh *bvh, uint32_t id,
                       const geomlib_trimesh *mesh)
{
    geomlib_bvh_node *n = &bvh->nodes[id];
    if (n->is_leaf) {
        uint32_t i;
        geomlib_aabb box = geomlib_aabb_empty();
        for (i = 0; i < n->data.leaf.primitive_count; ++i) {
            uint32_t tri = bvh->primitive_indices[n->data.leaf.primitive_offset+i];
            uint32_t a, b, c;
            geomlib_vec3 p;
            if (tri >= mesh->triangle_count) continue;
            a = mesh->indices[tri*3];
            b = mesh->indices[tri*3+1];
            c = mesh->indices[tri*3+2];
            if (a < mesh->vertex_count) {
                p = mesh->vertices[a];
                box = geomlib_aabb_expand(box, p);
            }
            if (b < mesh->vertex_count) {
                p = mesh->vertices[b];
                box = geomlib_aabb_expand(box, p);
            }
            if (c < mesh->vertex_count) {
                p = mesh->vertices[c];
                box = geomlib_aabb_expand(box, p);
            }
        }
        n->bounds = box;
        return;
    }
    refit_node(bvh, n->data.interior.left_child, mesh);
    refit_node(bvh, n->data.interior.right_child, mesh);
    n->bounds = geomlib_aabb_merge(bvh->nodes[n->data.interior.left_child].bounds,
                                   bvh->nodes[n->data.interior.right_child].bounds);
}

void geomlib_bvh_refit(geomlib_bvh *bvh, const geomlib_trimesh *mesh)
{
    if (!bvh || !mesh || !bvh->nodes || !mesh->vertices ||
        bvh->node_count == 0) return;
    bvh->mesh = mesh;
    refit_node(bvh, 0, mesh);
}
