/*
 * geomlib - Shape Diameter Function (see geomlib_shape_diameter.h)
 *
 * Clean-room implementation of Shlafman, Tal, Katz, "Symmetry of
 * Metallic Car Body Parts" / Shapira et al., "Shape Diameter Function", 2008:
 * cast rays in a cone around each vertex normal and take a robust statistic of
 * the hit distances, normalized by the model diagonal.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_shape_diameter.h"

#include <stdbool.h>
#include <math.h>

/* Double-sided Möller-Trumbore: the BVH ray tracer culls back faces, which
 * would make inward thickness rays miss every triangle of a closed surface. */
static bool hit_triangle(geomlib_vec3 org, geomlib_vec3 dir,
                         geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c,
                         float *t_out)
{
    geomlib_vec3 e1 = geomlib_vec3_sub(b,a);
    geomlib_vec3 e2 = geomlib_vec3_sub(c,a);
    geomlib_vec3 p = geomlib_vec3_cross(dir, e2);
    float det = geomlib_vec3_dot(e1, p);
    if (fabsf(det) < 1e-9f) return false;
    float inv = 1.0f/det;
    geomlib_vec3 tv = geomlib_vec3_sub(org, a);
    float u = geomlib_vec3_dot(tv, p)*inv;
    if (u < 0.0f || u > 1.0f) return false;
    geomlib_vec3 q = geomlib_vec3_cross(tv, e1);
    float v = geomlib_vec3_dot(dir, q)*inv;
    if (v < 0.0f || u+v > 1.0f) return false;
    float t = geomlib_vec3_dot(e2, q)*inv;
    if (t <= 0.0f) return false;
    *t_out = t;
    return true;
}

#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI_F 3.14159265358979323846f

geomlib_shape_diameter geomlib_shape_diameter_create(size_t count)
{
    geomlib_shape_diameter d;
    d.values = (float *)malloc((count?count:1)*sizeof(float));
    d.count = count;
    return d;
}

void geomlib_shape_diameter_free(geomlib_shape_diameter *d)
{
    if (!d) return;
    free(d->values);
    d->values = NULL;
    d->count = 0;
}

static int cmp_float(const void *a, const void *b)
{
    float x=*(const float*)a, y=*(const float*)b;
    return x<y ? -1 : (x>y ? 1 : 0);
}

int geomlib_compute_shape_diameter(const geomlib_trimesh *mesh,
                                   int ray_count,
                                   geomlib_shape_diameter *out)
{
    geomlib_aabb box;
    geomlib_vec3 size, diag_vec;
    float diag, eps;
    size_t V, i, k;
    int rays;
    float *samples;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (!out->values || out->count < V) {
        geomlib_shape_diameter_free(out);
        *out = geomlib_shape_diameter_create(V);
        if (!out->values) return GEOMLIB_MESH_ERROR_MEMORY;
    }
    rays = ray_count > 0 ? ray_count : 16;
    if (rays < 4) rays = 4;
    if (rays > 128) rays = 128;
    box = geomlib_trimesh_compute_aabb(mesh);
    size = geomlib_aabb_size(box);
    diag_vec = geomlib_vec3_create(size.x, size.y, size.z);
    diag = geomlib_vec3_length(diag_vec);
    if (!(diag > 1e-9f)) return GEOMLIB_MESH_ERROR_INVALID;
    eps = diag * 1e-4f;
    samples = (float *)malloc((size_t)rays*sizeof(float));
    if (!samples) return GEOMLIB_MESH_ERROR_MEMORY;
    geomlib_trimesh_compute_normals((geomlib_trimesh *)mesh);
    {
        geomlib_vec3 center = geomlib_aabb_center(box);
        for (i = 0; i < V; ++i) {
            if (geomlib_vec3_length(mesh->normals[i]) > 0.5f) continue;
            geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[i], center);
            if (geomlib_vec3_length(d) > 1e-9f)
                mesh->normals[i] = geomlib_vec3_normalize(d);
        }
    }

    for (i = 0; i < V; ++i) {
        geomlib_vec3 n = geomlib_vec3_normalize(mesh->normals[i]);
        geomlib_vec3 t1 = geomlib_vec3_sub(geomlib_vec3_create(1,0,0),
            geomlib_vec3_mul(n, geomlib_vec3_dot(geomlib_vec3_create(1,0,0), n)));
        geomlib_vec3 t2;
        size_t nhits = 0;
        if (geomlib_vec3_length(t1) < 1e-4f)
            t1 = geomlib_vec3_sub(geomlib_vec3_create(0,1,0),
                geomlib_vec3_mul(n, geomlib_vec3_dot(geomlib_vec3_create(0,1,0), n)));
        t1 = geomlib_vec3_normalize(t1);
        t2 = geomlib_vec3_cross(n, t1);
        for (k = 0; k < (size_t)rays; ++k) {
            float phi = 2.0f*PI_F*(float)k/(float)rays;
            float polar = 0.35f;
            /* Rays point into the solid (along -normal) to measure local
             * thickness; the offset origin keeps them clear of the source
             * triangle. */
            geomlib_vec3 dir = geomlib_vec3_add(
                geomlib_vec3_mul(n, -cosf(polar)),
                geomlib_vec3_add(geomlib_vec3_mul(t1, sinf(polar)*cosf(phi)),
                                 geomlib_vec3_mul(t2, sinf(polar)*sinf(phi))));
            /* Start slightly inside the surface so the source triangle is
             * behind the ray; add the offset back to the measured distance. */
            geomlib_vec3 origin = geomlib_vec3_sub(mesh->vertices[i],
                                                   geomlib_vec3_mul(n, diag*0.005f));
            {
                geomlib_vec3 d = geomlib_vec3_normalize(dir);
                float best = FLT_MAX;
                size_t tri;
                for (tri = 0; tri < mesh->triangle_count; ++tri) {
                    float t;
                    uint32_t i0=mesh->indices[tri*3], i1=mesh->indices[tri*3+1],
                             i2=mesh->indices[tri*3+2];
                    if (!hit_triangle(origin, d, mesh->vertices[i0],
                                      mesh->vertices[i1], mesh->vertices[i2], &t)) continue;
                    if (t > eps && t < best) best = t;
                }
                if (best < FLT_MAX) samples[nhits++] = best + diag*0.005f;
            }
        }
        if (nhits) {
            qsort(samples, nhits, sizeof(float), cmp_float);
            out->values[i] = samples[nhits/2] / diag;
        } else {
            out->values[i] = 0.0f;
        }
    }
    free(samples);
    return GEOMLIB_MESH_SUCCESS;
}
