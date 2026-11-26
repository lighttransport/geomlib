/*
 * geomlib - Triangle-mesh intersection (see geomlib_intersection.h)
 *
 * Clean-room implementation. An intersection edge lies at the crossing of one
 * face of A and one face of B. Every vertex of the intersection polyline is
 * either an A-edge/B-face crossing or a B-edge/A-face crossing, so we trace
 * edges of both meshes through a BVH of the other, tag each crossing with its
 * (A face, B face) pair, sort angularly inside each pair, and join the loops.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_intersection.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_intersection geomlib_intersection_create(void)
{
    geomlib_intersection r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_intersection_free(geomlib_intersection *r)
{
    if (!r) return;
    free(r->segments);
    memset(r, 0, sizeof(*r));
}

typedef struct { uint32_t aface, bface; geomlib_vec3 p; } crossing;

static int cmp_crossing(const void *pa, const void *pb)
{
    const crossing *a=(const crossing*)pa, *b=(const crossing*)pb;
    if (a->aface < b->aface) return -1;
    if (a->aface > b->aface) return 1;
    if (a->bface < b->bface) return -1;
    if (a->bface > b->bface) return 1;
    return 0;
}

static int close_point(geomlib_vec3 a, geomlib_vec3 b)
{
    return geomlib_vec3_length(geomlib_vec3_sub(a,b)) < 1e-4f;
}

/* Which faces of `m` contain both endpoints of edge (u,v)? */
static size_t incident_faces(const geomlib_trimesh *m, uint32_t u, uint32_t v,
                             uint32_t *faces, size_t max)
{
    size_t n = 0, t;
    for (t = 0; t < m->triangle_count && n < max; ++t) {
        uint32_t x=m->indices[t*3], y=m->indices[t*3+1], z=m->indices[t*3+2];
        if (!((x==u&&y==v)||(x==v&&y==u)||(x==u&&z==v)||
              (x==v&&z==u)||(y==u&&z==v)||(y==v&&z==u))) continue;
        faces[n++] = (uint32_t)t;
    }
    return n;
}

static int add_crossing(crossing **arr, size_t *n, size_t *cap,
                        uint32_t aface, uint32_t bface, geomlib_vec3 p)
{
    size_t k;
    for (k = *n; k > 0 && (*arr)[k-1].aface == aface &&
                       (*arr)[k-1].bface == bface; --k)
        if (close_point((*arr)[k-1].p, p)) return 0;
    if (*n == *cap) {
        size_t nc = *cap ? (*cap)*2 : 128;
        crossing *na = (crossing *)realloc(*arr, nc*sizeof(crossing));
        if (!na) return -1;
        *arr = na; *cap = nc;
    }
    (*arr)[*n].aface = aface;
    (*arr)[*n].bface = bface;
    (*arr)[*n].p = p;
    (*n)++;
    return 0;
}

/* Trace every edge of `src` against BVH of `dst`, tagging with src faces. */
static int inside_box(geomlib_vec3 p, geomlib_aabb b, float eps)
{
    return p.x >= b.min.x-eps && p.x <= b.max.x+eps &&
           p.y >= b.min.y-eps && p.y <= b.max.y+eps &&
           p.z >= b.min.z-eps && p.z <= b.max.z+eps;
}

static int trace_edges(const geomlib_trimesh *src, const geomlib_bvh *dst_bvh,
                       int src_is_a, crossing **cross, size_t *ncross,
                       size_t *ccap, geomlib_aabb boxa, geomlib_aabb boxb)
{
    size_t t;
    for (t = 0; t < src->triangle_count; ++t) {
        uint32_t f[3];
        size_t k;
        for (k = 0; k < 3; ++k) f[k] = src->indices[t*3+k];
        for (k = 0; k < 3; ++k) {
            uint32_t u = f[k], v = f[(k+1)%3];
            uint32_t faces[8];
            size_t nf, fi, iter;
            geomlib_vec3 p, q, delta, dir, o;
            float len;
            if (u >= src->vertex_count || v >= src->vertex_count) continue;
            p = src->vertices[u]; q = src->vertices[v];
            nf = incident_faces(src, u, v, faces, 8);
            if (!nf) continue;
            delta = geomlib_vec3_sub(q,p);
            len = geomlib_vec3_length(delta);
            if (!(len > 1e-6f)) continue;
            dir = geomlib_vec3_mul(delta, 1.0f/len);
            o = p;
            for (iter = 0; iter < 32; ++iter) {
                geomlib_ray ray = geomlib_ray_create(o, dir);
                geomlib_ray_hit hit;
                ray.t_min = 1e-5f;
                ray.t_max = len - 1e-6f;
                if (ray.t_max <= ray.t_min ||
                    !geomlib_bvh_trace_ray(dst_bvh, &ray, &hit)) break;
                /* A true intersection point lies within both AABBs. */
                if (!inside_box(hit.position, boxa, 1e-4f) ||
                    !inside_box(hit.position, boxb, 1e-4f)) continue;
                for (fi = 0; fi < nf; ++fi) {
                    uint32_t af = src_is_a ? faces[fi] : hit.primitive_id;
                    uint32_t bf = src_is_a ? hit.primitive_id : faces[fi];
                    if (add_crossing(cross, ncross, ccap, af, bf, hit.position) < 0)
                        return -1;
                }
                o = geomlib_vec3_add(hit.position, geomlib_vec3_mul(dir, 1e-4f));
            }
        }
    }
    return 0;
}

typedef struct { crossing c; double angle; } angular;
static int cmp_angular(const void *pa, const void *pb)
{
    const angular *a=(const angular*)pa, *b=(const angular*)pb;
    return a->angle < b->angle ? -1 : (a->angle > b->angle ? 1 : 0);
}

static int add_segment(geomlib_intersection *r, geomlib_vec3 a, geomlib_vec3 b,
                       size_t *cap)
{
    if (close_point(a,b)) return 0;
    if (r->segment_count == *cap) {
        size_t nc = *cap ? (*cap)*2 : 64;
        geomlib_vec3 *ns = (geomlib_vec3 *)realloc(r->segments,
                                                   nc*2*sizeof(geomlib_vec3));
        if (!ns) return -1;
        r->segments = ns; *cap = nc;
    }
    r->segments[r->segment_count*2] = a;
    r->segments[r->segment_count*2+1] = b;
    r->segment_count++;
    return 0;
}

geomlib_intersection geomlib_mesh_mesh_intersection(const geomlib_trimesh *a,
                                                     const geomlib_trimesh *b)
{
    geomlib_intersection result = geomlib_intersection_create();
    geomlib_bvh bvha, bvhb;
    crossing *cross = NULL;
    angular *ring = NULL;
    size_t ncross = 0, ccap = 0, cap = 0, i;

    if (!a || !b || !a->vertices || !a->indices ||
        !b->vertices || !b->indices || !a->triangle_count || !b->triangle_count)
        return result;
    if (geomlib_bvh_build(&bvha, a, GEOMLIB_BVH_QUALITY_MEDIUM) !=
        GEOMLIB_MESH_SUCCESS) return result;
    if (geomlib_bvh_build(&bvhb, b, GEOMLIB_BVH_QUALITY_MEDIUM) !=
        GEOMLIB_MESH_SUCCESS) { geomlib_bvh_free(&bvha); return result; }

    {
        geomlib_aabb boxa = geomlib_trimesh_compute_aabb(a);
        geomlib_aabb boxb = geomlib_trimesh_compute_aabb(b);
        if (trace_edges(a, &bvhb, 1, &cross, &ncross, &ccap, boxa, boxb) < 0 ||
            trace_edges(b, &bvha, 0, &cross, &ncross, &ccap, boxa, boxb) < 0)
            goto cleanup;
    }
    geomlib_bvh_free(&bvha);
    geomlib_bvh_free(&bvhb);
    if (!ncross) goto cleanup;

    qsort(cross, ncross, sizeof(crossing), cmp_crossing);
    ring = (angular *)malloc(ncross*sizeof(angular));
    if (!ring) goto cleanup;

    for (i = 0; i < ncross;) {
        size_t j = i+1, m, n;
        geomlib_vec3 center = geomlib_vec3_create(0,0,0);
        geomlib_vec3 nrm = geomlib_vec3_create(0,0,0), t1, t2, helper;
        while (j < ncross && cross[j].aface == cross[i].aface &&
                             cross[j].bface == cross[i].bface) ++j;
        for (m = i; m < j; ++m) center = geomlib_vec3_add(center, cross[m].p);
        center = geomlib_vec3_mul(center, 1.0f/(float)(j-i));
        {
            uint32_t f = cross[i].aface;
            if (f < a->triangle_count) {
                uint32_t i0=a->indices[f*3], i1=a->indices[f*3+1], i2=a->indices[f*3+2];
                if (i0<a->vertex_count && i1<a->vertex_count && i2<a->vertex_count) {
                    geomlib_vec3 e1=geomlib_vec3_sub(a->vertices[i1],a->vertices[i0]);
                    geomlib_vec3 e2=geomlib_vec3_sub(a->vertices[i2],a->vertices[i0]);
                    nrm = geomlib_vec3_cross(e1,e2);
                }
            }
        }
        if (geomlib_vec3_length(nrm) < 1e-12f) nrm = geomlib_vec3_create(0,0,1);
        nrm = geomlib_vec3_normalize(nrm);
        helper = fabsf(nrm.x) < 0.9f ? geomlib_vec3_create(1,0,0)
                                    : geomlib_vec3_create(0,1,0);
        t1 = geomlib_vec3_normalize(geomlib_vec3_cross(nrm, helper));
        t2 = geomlib_vec3_cross(nrm, t1);
        for (m = i; m < j; ++m) {
            geomlib_vec3 d = geomlib_vec3_sub(cross[m].p, center);
            ring[m-i].c = cross[m];
            ring[m-i].angle = atan2(geomlib_vec3_dot(d,t2), geomlib_vec3_dot(d,t1));
        }
        qsort(ring, j-i, sizeof(angular), cmp_angular);
        n = j-i;
        if (n == 2) {
            if (add_segment(&result, ring[0].c.p, ring[1].c.p, &cap) < 0) goto cleanup;
        } else {
            for (m = 0; m < n; ++m) {
                size_t nx = (m+1)%n;
                if (add_segment(&result, ring[m].c.p, ring[nx].c.p, &cap) < 0) goto cleanup;
            }
        }
        i = j;
    }

cleanup:
    free(cross);
    free(ring);
    return result;
}
