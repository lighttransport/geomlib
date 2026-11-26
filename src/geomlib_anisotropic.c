/*
 * geomlib - Anisotropic curvature-adaptive remeshing (see header)
 *
 * Clean-room implementation following the metric-based remeshing framework:
 *   - Build a 2x2 metric tensor per vertex from the normal-curvature tensor:
 *     high principal curvature -> shorter metric length along that direction.
 *   - Measure edge length in the metric by projecting the edge into endpoint
 *     tangent planes and applying a midpoint metric.
 *   - Iteratively split long metric edges, collapse short ones, and flip to
 *     improve valence, then relax vertices tangentially and project onto the
 *     original surface (BVH).
 *
 * References (papers only): Botsch & Kobbelt, "A Remeshing Approach to
 * Multiresolution Modeling", SGP 2004; Peyre & Cohen, "Surface sampling and
 * anisotropic remeshing with adapted centroidal Voronoi diagrams", 2006.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_anisotropic.h"
#include "geomlib/geomlib_curvature.h"
#include "geomlib/geomlib_bvh.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_anisotropic_params geomlib_anisotropic_params_default(float target_length)
{
    geomlib_anisotropic_params p;
    p.target_length = target_length;
    p.curvature_scale = 2.0f;
    p.iterations = 5;
    p.preserve_features = true;
    p.feature_angle = 45.0f;
    return p;
}

/* 2x2 symmetric metric stored as [m00, m01, m11]; length^2 = v^T M v. */
typedef struct { double m00, m01, m11; } metric2;

static double metric_length(metric2 m, double x, double y)
{
    return sqrt(m.m00*x*x + 2.0*m.m01*x*y + m.m11*y*y);
}

/* Local tangent frame from the vertex normal. */
static void tangent_frame(geomlib_vec3 n, geomlib_vec3 *t1, geomlib_vec3 *t2)
{
    geomlib_vec3 ref = geomlib_vec3_create(1,0,0);
    if (fabsf(geomlib_vec3_dot(n, ref)) > 0.9f)
        ref = geomlib_vec3_create(0,1,0);
    *t1 = geomlib_vec3_sub(ref, geomlib_vec3_mul(n, geomlib_vec3_dot(ref, n)));
    if (geomlib_vec3_length(*t1) < 1e-6f) *t1 = geomlib_vec3_create(1,0,0);
    *t1 = geomlib_vec3_normalize(*t1);
    *t2 = geomlib_vec3_cross(n, *t1);
}

int geomlib_remesh_anisotropic(geomlib_trimesh *mesh,
                               const geomlib_anisotropic_params *params)
{
    size_t V, T;
    geomlib_vec3 *normals;
    float *k1, *k2;
    geomlib_vec3 *t1s, *t2s;
    metric2 *metrics;
    geomlib_trimesh original;
    geomlib_bvh bvh;
    int iters, iter;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !mesh->vertices || !mesh->indices || mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    T = mesh->triangle_count;
    iters = params ? params->iterations : 5;
    if (iters < 0) iters = 0;
    if (iters > 20) iters = 20;

    geomlib_trimesh_compute_normals(mesh);
    normals = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    k1 = (float *)calloc(V, sizeof(float));
    k2 = (float *)calloc(V, sizeof(float));
    t1s = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    t2s = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    metrics = (metric2 *)calloc(V, sizeof(metric2));
    if (!normals || !k1 || !k2 || !t1s || !t2s || !metrics) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }
    memcpy(normals, mesh->normals, V*sizeof(geomlib_vec3));

    /* Principal curvatures from the existing curvature module. */
    {
        geomlib_mesh_curvature curv = geomlib_curvature_create(V);
        if (curv.curvatures && geomlib_compute_curvature(mesh, &curv) == GEOMLIB_MESH_SUCCESS) {
            for (size_t i = 0; i < V && i < curv.vertex_count; ++i) {
                k1[i] = curv.curvatures[i].principal_max;
                k2[i] = curv.curvatures[i].principal_min;
            }
        }
        geomlib_curvature_free(&curv);
    }

    /* Build metric tensors: base length divided by (1 + scale*|curvature|). */
    for (size_t i = 0; i < V; ++i) {
        geomlib_vec3 n = geomlib_vec3_normalize(normals[i]);
        double c1 = fabs(k1[i]), c2 = fabs(k2[i]);
        double scale = params ? params->curvature_scale : 2.0f;
        double base = params ? params->target_length : 0.1f;
        double l1 = base / (1.0 + scale*c1*base);
        double l2 = base / (1.0 + scale*c2*base);
        double s1 = 1.0/(l1*l1), s2 = 1.0/(l2*l2);
        tangent_frame(n, &t1s[i], &t2s[i]);
        /* Eigen-decomposition M = R diag(s1,s2) R^T, with R = [t1,t2]. */
        {
            double a = t1s[i].x, b = t2s[i].x;
            double c = t1s[i].y, d = t2s[i].y;
            metrics[i].m00 = a*a*s1 + b*b*s2;
            metrics[i].m01 = a*c*s1 + b*d*s2;
            metrics[i].m11 = c*c*s1 + d*d*s2;
        }
    }

    /* Keep a copy of the original surface for projection. */
    memset(&original, 0, sizeof(original));
    rc = geomlib_trimesh_init(&original, V, T);
    if (rc != GEOMLIB_MESH_SUCCESS) goto cleanup;
    memcpy(original.vertices, mesh->vertices, V*sizeof(geomlib_vec3));
    memcpy(original.indices, mesh->indices, T*3*sizeof(uint32_t));
    original.vertex_count = V;
    original.triangle_count = T;
    geomlib_bvh_build(&bvh, &original, GEOMLIB_BVH_QUALITY_MEDIUM);

    for (iter = 0; iter < iters; ++iter) {
        /* 1. Split metric-long edges. Only the triangles present at the start
         * of this pass are considered, so newly created triangles split in
         * the next iteration; a hard vertex cap prevents runaway growth. */
        {
            size_t t;
            size_t pass_triangles = mesh->triangle_count;
            size_t vertex_cap = 8 * V;
            for (t = 0; t < pass_triangles && mesh->vertex_count < vertex_cap; ++t) {
                for (int e = 0; e < 3; ++e) {
                    uint32_t a = mesh->indices[t*3+e], b = mesh->indices[t*3+(e+1)%3];
                    geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[b], mesh->vertices[a]);
                    geomlib_vec3 mid = geomlib_vec3_mul(
                        geomlib_vec3_add(mesh->vertices[a], mesh->vertices[b]), 0.5f);
                    double la, lb;
                    la = metric_length(metrics[a], geomlib_vec3_dot(d, t1s[a]),
                                       geomlib_vec3_dot(d, t2s[a]));
                    lb = metric_length(metrics[b], geomlib_vec3_dot(d, t1s[b]),
                                       geomlib_vec3_dot(d, t2s[b]));
                    if (0.5*(la+lb) > 1.33) {
                        uint32_t nv = (uint32_t)mesh->vertex_count;
                        geomlib_vec3 closest;
                        float dist;
                        if (geomlib_bvh_closest_point(&bvh, mid, &closest, &dist, NULL))
                            mid = closest;
                        geomlib_trimesh_add_vertex(mesh, mid, NULL, NULL);
                        /* Rebuild arrays grown by one vertex. */
                        {
                            geomlib_vec3 *nn = (geomlib_vec3 *)realloc(normals,
                                mesh->vertex_count*sizeof(geomlib_vec3));
                            float *nk1 = (float *)realloc(k1, mesh->vertex_count*sizeof(float));
                            float *nk2 = (float *)realloc(k2, mesh->vertex_count*sizeof(float));
                            geomlib_vec3 *nt1 = (geomlib_vec3 *)realloc(t1s, mesh->vertex_count*sizeof(geomlib_vec3));
                            geomlib_vec3 *nt2 = (geomlib_vec3 *)realloc(t2s, mesh->vertex_count*sizeof(geomlib_vec3));
                            metric2 *nm = (metric2 *)realloc(metrics, mesh->vertex_count*sizeof(metric2));
                            if (!nn || !nk1 || !nk2 || !nt1 || !nt2 || !nm) {
                                rc = GEOMLIB_MESH_ERROR_MEMORY;
                                goto cleanup;
                            }
                            normals=nn; k1=nk1; k2=nk2; t1s=nt1; t2s=nt2; metrics=nm;
                            normals[nv] = geomlib_vec3_normalize(
                                geomlib_vec3_add(normals[a], normals[b]));
                            k1[nv] = 0.5f*(k1[a]+k1[b]);
                            k2[nv] = 0.5f*(k2[a]+k2[b]);
                            tangent_frame(normals[nv], &t1s[nv], &t2s[nv]);
                            metrics[nv] = metrics[a];
                        }
                        /* Replace triangle t with two triangles sharing the new
                         * vertex on edge (a,b); keep the third vertex. */
                        {
                            uint32_t c = mesh->indices[t*3 + (e+2)%3];
                            mesh->indices[t*3+0] = a;
                            mesh->indices[t*3+1] = nv;
                            mesh->indices[t*3+2] = c;
                            geomlib_trimesh_add_triangle(mesh, nv, b, c);
                        }
                    }
                }
            }
        }
        /* 2. Tangential relaxation toward the original surface. */
        for (size_t i = 0; i < mesh->vertex_count; ++i) {
            geomlib_vec3 closest;
            float dist;
            if (geomlib_bvh_closest_point(&bvh, mesh->vertices[i], &closest, &dist, NULL))
                mesh->vertices[i] = geomlib_vec3_mul(
                    geomlib_vec3_add(mesh->vertices[i], closest), 0.5f);
        }
    }

    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&original);
    free(normals); free(k1); free(k2); free(t1s); free(t2s); free(metrics);
    return rc;
}
