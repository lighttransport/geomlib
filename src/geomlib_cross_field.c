/*
 * geomlib - Cross fields (see geomlib_cross_field.h)
 *
 * Clean-room implementation of the standard cross-field combing procedure:
 * seed directions from principal curvature, propagate across mesh edges with
 * the 90-degree symmetry quotient taken into account.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_cross_field.h"
#include "geomlib/geomlib_curvature.h"
#include "geomlib/geomlib_topology.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI_F 3.14159265358979323846f

geomlib_cross_field geomlib_cross_field_create(size_t count)
{
    geomlib_cross_field f;
    f.angles = (float *)malloc((count?count:1)*sizeof(float));
    f.normals_x = (float *)malloc((count?count:1)*sizeof(float));
    f.normals_y = (float *)malloc((count?count:1)*sizeof(float));
    f.count = count;
    if (!f.angles || !f.normals_x || !f.normals_y) geomlib_cross_field_free(&f);
    return f;
}

void geomlib_cross_field_free(geomlib_cross_field *f)
{
    if (!f) return;
    free(f->angles);
    free(f->normals_x);
    free(f->normals_y);
    f->angles = NULL;
    f->normals_x = NULL;
    f->normals_y = NULL;
    f->count = 0;
}

static float wrap90(float a)
{
    float q = PI_F*0.5f;
    while (a >= q) a -= q;
    while (a < 0.0f) a += q;
    return a;
}

int geomlib_compute_cross_field(const geomlib_trimesh *mesh,
                                geomlib_cross_field *out)
{
    geomlib_vertex_adjacency adj = {0};
    geomlib_mesh_curvature curv;
    size_t V, i, j, iter;
    geomlib_vec3 *nrm, *t1s, *t2s;
    bool *seen;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    if (!out->angles || out->count < V) {
        geomlib_cross_field_free(out);
        *out = geomlib_cross_field_create(V);
        if (!out->angles) return GEOMLIB_MESH_ERROR_MEMORY;
    }

    geomlib_trimesh_compute_normals((geomlib_trimesh *)mesh);
    nrm = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    t1s = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    t2s = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    seen = (bool *)calloc(V, sizeof(bool));
    if (!nrm || !t1s || !t2s || !seen) {
        free(nrm); free(t1s); free(t2s); free(seen);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < V; ++i) nrm[i] = geomlib_vec3_normalize(mesh->normals[i]);

    /* Stable tangent frame per vertex. */
    for (i = 0; i < V; ++i) {
        geomlib_vec3 ref = geomlib_vec3_create(1,0,0);
        geomlib_vec3 t1, t2;
        if (fabsf(geomlib_vec3_dot(nrm[i], ref)) > 0.9f)
            ref = geomlib_vec3_create(0,1,0);
        t1 = geomlib_vec3_sub(ref, geomlib_vec3_mul(nrm[i],
                              geomlib_vec3_dot(ref, nrm[i])));
        if (geomlib_vec3_length(t1) < 1e-6f) t1 = geomlib_vec3_create(1,0,0);
        t1 = geomlib_vec3_normalize(t1);
        t2 = geomlib_vec3_cross(nrm[i], t1);
        t1s[i] = t1;
        t2s[i] = t2;
        out->normals_x[i] = t1.x;
        out->normals_y[i] = t2.x;
        out->angles[i] = 0.0f;
    }

    /* Seed from principal curvature directions. */
    curv = geomlib_curvature_create(V);
    if (curv.curvatures && geomlib_compute_curvature(mesh, &curv) == GEOMLIB_MESH_SUCCESS) {
        for (i = 0; i < V && i < curv.vertex_count; ++i) {
            geomlib_vec3 d = curv.curvatures[i].principal_dir_max;
            float u = geomlib_vec3_dot(d, t1s[i]);
            float v = geomlib_vec3_dot(d, t2s[i]);
            out->angles[i] = wrap90(atan2f(v, u));
        }
    }
    geomlib_curvature_free(&curv);

    /* Comb by propagating angles across edges with 90-degree symmetry. */
    adj = geomlib_build_vertex_adjacency(mesh);
    if (!adj.neighbor_counts) {
        free(nrm); free(t1s); free(t2s); free(seen);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    (void)seen;
    /* Relax every edge: replace each angle by the congruent representative
     * (mod 90 degrees) nearest to its neighbour's angle. */
    for (iter = 0; iter < V; ++iter) {
        int changed = 0;
        for (i = 0; i < V; ++i) {
            for (j = 0; j < adj.neighbor_counts[i]; ++j) {
                uint32_t nb = adj.neighbors[i][j];
                float q = PI_F*0.5f;
                float a = out->angles[nb];
                float target = out->angles[i];
                float k = floorf((target - a)/q + 0.5f);
                float cand = a + k*q;
                if (fabsf(cand-a) > 1e-6f) {
                    out->angles[nb] = cand;
                    changed = 1;
                }
            }
        }
        if (!changed) break;
    }

    geomlib_vertex_adjacency_free(&adj);
    free(nrm); free(t1s); free(t2s);
    free(seen);
    return GEOMLIB_MESH_SUCCESS;
}
