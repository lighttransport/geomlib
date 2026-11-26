/*
 * geomlib - Vertex-clustering decimation (see geomlib_clustering.h)
 *
 * Clean-room implementation of uniform vertex clustering: hash vertices by
 * voxel cell, average each cluster's positions, remap triangle indices to
 * cluster ids and drop degenerate/duplicate triangles.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_clustering.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int x, y, z;
    uint32_t id;
    double sx, sy, sz;
    size_t count;
} cluster;

static int cmp_cluster(const void *pa, const void *pb)
{
    const cluster *a = (const cluster *)pa, *b = (const cluster *)pb;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    if (a->z < b->z) return -1;
    if (a->z > b->z) return 1;
    return 0;
}

typedef struct {
    uint32_t a, b, c;      /* original winding */
    uint32_t s0, s1, s2;   /* sorted key for duplicate detection */
} tkey;

static int cmp_tri(const void *pa, const void *pb)
{
    const tkey *x = (const tkey *)pa, *y = (const tkey *)pb;
    if (x->s0 < y->s0) return -1;
    if (x->s0 > y->s0) return 1;
    if (x->s1 < y->s1) return -1;
    if (x->s1 > y->s1) return 1;
    if (x->s2 < y->s2) return -1;
    if (x->s2 > y->s2) return 1;
    return 0;
}

int geomlib_decimate_vertex_clustering(const geomlib_trimesh *input,
                                       float voxel_size,
                                       geomlib_trimesh *output)
{
    cluster *clusters;
    uint32_t *remap;
    geomlib_vec3 *positions;
    size_t n, i, count, ntri, out_tri;
    tkey *tris;

    if (!input || !output || !input->vertices || !input->indices ||
        input->triangle_count == 0 || !(voxel_size > 0.0f))
        return GEOMLIB_MESH_ERROR_INVALID;
    n = input->vertex_count;
    clusters = (cluster *)malloc(n*sizeof(cluster));
    remap = (uint32_t *)malloc(n*sizeof(uint32_t));
    positions = (geomlib_vec3 *)malloc(n*sizeof(geomlib_vec3));
    tris = (tkey *)malloc(input->triangle_count*sizeof(tkey));
    if (!clusters || !remap || !positions || !tris) {
        free(clusters); free(remap); free(positions); free(tris);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (i = 0; i < n; ++i) {
        geomlib_vec3 p = input->vertices[i];
        clusters[i].x = (int)floorf(p.x/voxel_size);
        clusters[i].y = (int)floorf(p.y/voxel_size);
        clusters[i].z = (int)floorf(p.z/voxel_size);
        clusters[i].id = 0;
        clusters[i].sx = p.x;
        clusters[i].sy = p.y;
        clusters[i].sz = p.z;
        clusters[i].count = 1;
    }
    qsort(clusters, n, sizeof(cluster), cmp_cluster);
    count = 0;
    for (i = 0; i < n;) {
        size_t j = i+1;
        double sx = clusters[i].sx, sy = clusters[i].sy, sz = clusters[i].sz;
        size_t c = clusters[i].count;
        while (j < n && cmp_cluster(&clusters[i], &clusters[j]) == 0) {
            sx += clusters[j].sx;
            sy += clusters[j].sy;
            sz += clusters[j].sz;
            c += clusters[j].count;
            ++j;
        }
        clusters[i].sx = sx;
        clusters[i].sy = sy;
        clusters[i].sz = sz;
        clusters[i].count = c;
        clusters[i].id = (uint32_t)count++;
        i = j;
    }

    /* Restore original vertex order in the remap table by matching cells. */
    for (i = 0; i < n; ++i) {
        geomlib_vec3 p = input->vertices[i];
        cluster key;
        cluster *found;
        key.x = (int)floorf(p.x/voxel_size);
        key.y = (int)floorf(p.y/voxel_size);
        key.z = (int)floorf(p.z/voxel_size);
        key.id = 0; key.sx=key.sy=key.sz=0; key.count=0;
        found = (cluster *)bsearch(&key, clusters, n, sizeof(cluster), cmp_cluster);
        remap[i] = found ? found->id : 0;
    }
    for (i = 0; i < count; ++i) {
        positions[i] = geomlib_vec3_create(
            (float)(clusters[i].sx/clusters[i].count),
            (float)(clusters[i].sy/clusters[i].count),
            (float)(clusters[i].sz/clusters[i].count));
    }

    ntri = 0;
    for (i = 0; i < input->triangle_count; ++i) {
        uint32_t a = remap[input->indices[i*3]];
        uint32_t b = remap[input->indices[i*3+1]];
        uint32_t c = remap[input->indices[i*3+2]];
        if (a == b || b == c || a == c) continue;
        tkey key;
        uint32_t s0=a, s1=b, s2=c, tmp;
        if (s0 > s1) { tmp=s0; s0=s1; s1=tmp; }
        if (s1 > s2) { tmp=s1; s1=s2; s2=tmp; }
        if (s0 > s1) { tmp=s0; s0=s1; s1=tmp; }
        key.a = remap[input->indices[i*3]];
        key.b = remap[input->indices[i*3+1]];
        key.c = remap[input->indices[i*3+2]];
        key.s0 = s0; key.s1 = s1; key.s2 = s2;
        tris[ntri++] = key;
    }
    qsort(tris, ntri, sizeof(tkey), cmp_tri);
    out_tri = 0;
    for (i = 0; i < ntri; ++i) {
        if (out_tri && cmp_tri(&tris[out_tri-1], &tris[i]) == 0) continue;
        tris[out_tri++] = tris[i];
    }

    {
        geomlib_trimesh result;
        int rc = geomlib_trimesh_init(&result, count, out_tri);
        if (rc != GEOMLIB_MESH_SUCCESS) {
            free(clusters); free(remap); free(positions); free(tris);
            return rc;
        }
        for (i = 0; i < count; ++i)
            geomlib_trimesh_add_vertex(&result, positions[i], NULL, NULL);
        for (i = 0; i < out_tri; ++i)
            geomlib_trimesh_add_triangle(&result, tris[i].a, tris[i].b, tris[i].c);
        *output = result;
    }
    free(clusters); free(remap); free(positions); free(tris);
    return GEOMLIB_MESH_SUCCESS;
}
