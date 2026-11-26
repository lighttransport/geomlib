/*
 * geomlib - ARAP surface deformation (see geomlib_deform.h)
 *
 * Clean-room implementation of Sorkine & Alexa, "As-Rigid-As-Possible Surface
 * Modeling", Eurographics Symposium on Geometry Processing 2007.
 *
 * Local step: best-fit rotation per vertex from rest edges to current edges
 * (Horn quaternion method). Global step: sparse least-squares solve for the
 * new positions with handle positions pinned. C11, wasm32 safe.
 */

#include "geomlib/geomlib_deform.h"
#include "geomlib/geomlib_sparse.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t a, b; } edge_key;

static int cmp_edge(const void *pa, const void *pb)
{
    const edge_key *x = (const edge_key *)pa;
    const edge_key *y = (const edge_key *)pb;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

/* Best-fit 3x3 rotation (Horn quaternion) aligning centered sets a -> b. */
static void fit_rotation(const geomlib_vec3 *a, const geomlib_vec3 *b,
                         size_t n, double out[9])
{
    double H[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double N[16] = {0}, eig[4], vec[16];
    size_t i;
    double sxx,syy,szz,sxy,syx,sxz,szx,syz,szy;

    for (i = 0; i < n; ++i) {
        double ax=a[i].x, ay=a[i].y, az=a[i].z;
        double bx=b[i].x, by=b[i].y, bz=b[i].z;
        H[0][0]+=ax*bx; H[0][1]+=ax*by; H[0][2]+=ax*bz;
        H[1][0]+=ay*bx; H[1][1]+=ay*by; H[1][2]+=ay*bz;
        H[2][0]+=az*bx; H[2][1]+=az*by; H[2][2]+=az*bz;
    }
    sxx=H[0][0]; syy=H[1][1]; szz=H[2][2];
    sxy=H[0][1]; syx=H[1][0]; sxz=H[0][2]; szx=H[2][0];
    syz=H[1][2]; szy=H[2][1];
    N[0]=sxx+syy+szz; N[5]=sxx-syy-szz; N[10]=-sxx+syy-szz; N[15]=-sxx-syy+szz;
    N[1]=N[4]=syz-szy; N[2]=N[8]=szx-sxz; N[3]=N[12]=sxy-syx;
    N[6]=N[9]=sxy+syx; N[7]=N[13]=sxz+szx; N[11]=N[14]=syz+szy;
    if (geomlib_eigen_sym_dense(N, 4, eig, vec) != GEOMLIB_SPARSE_SUCCESS) {
        out[0]=out[4]=out[8]=1.0;
        return;
    }
    {
        const double *v = vec + 12;
        double x=v[0], y=v[1], z=v[2], w=v[3];
        double norm = sqrt(w*w+x*x+y*y+z*z);
        if (norm < 1e-12) { out[0]=out[4]=out[8]=1.0; return; }
        x/=norm; y/=norm; z/=norm; w/=norm;
        out[0]=w*w+x*x-y*y-z*z; out[1]=2*(x*y-w*z); out[2]=2*(x*z+w*y);
        out[3]=2*(x*y+w*z); out[4]=w*w-x*x+y*y-z*z; out[5]=2*(y*z-w*x);
        out[6]=2*(x*z-w*y); out[7]=2*(y*z+w*x); out[8]=w*w-x*x-y*y+z*z;
    }
}

int geomlib_deform_arap(const geomlib_trimesh *input,
                        const geomlib_deform_handle *handles,
                        size_t handle_count,
                        int iterations,
                        geomlib_trimesh *output)
{
    size_t V, T, E = 0, i, k, iter;
    edge_key *edges;
    size_t *vcount, *voff;
    uint32_t *vneighbors;
    double *rot;
    double *curx, *cury, *curz;
    geomlib_vec3 *rest;
    uint32_t *rows = NULL, *cols = NULL;
    double *vals = NULL;
    double *solx = NULL, *soly = NULL, *solz = NULL;
    geomlib_sparse_solve_params sp;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!input || !output || !input->vertices || !input->indices ||
        input->triangle_count == 0 || !handles || handle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = input->vertex_count;
    T = input->triangle_count;
    if (iterations <= 0) iterations = 10;
    if (iterations > 100) iterations = 100;

    edges = (edge_key *)malloc(3*T*sizeof(edge_key));
    vcount = (size_t *)calloc(V, sizeof(size_t));
    if (!edges || !vcount) { free(edges); free(vcount); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < T; ++i) {
        for (k = 0; k < 3; ++k) {
            uint32_t a = input->indices[i*3+k], b2 = input->indices[i*3+(k+1)%3];
            edges[E].a = a < b2 ? a : b2;
            edges[E].b = a < b2 ? b2 : a;
            E++;
        }
    }
    qsort(edges, E, sizeof(edge_key), cmp_edge);
    {
        size_t n = 0;
        for (i = 0; i < E;) {
            size_t j = i+1;
            while (j < E && cmp_edge(&edges[i], &edges[j]) == 0) ++j;
            edges[n++] = edges[i];
            i = j;
        }
        E = n;
    }
    for (i = 0; i < E; ++i) { vcount[edges[i].a]++; vcount[edges[i].b]++; }
    voff = (size_t *)malloc((V+1)*sizeof(size_t));
    vneighbors = (uint32_t *)malloc(2*E*sizeof(uint32_t));
    rest = (geomlib_vec3 *)malloc(V*sizeof(geomlib_vec3));
    rot = (double *)calloc(V*9, sizeof(double));
    curx = (double *)calloc(V, sizeof(double));
    cury = (double *)calloc(V, sizeof(double));
    curz = (double *)calloc(V, sizeof(double));
    if (!voff || !vneighbors || !rest || !rot || !curx || !cury || !curz) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }
    voff[0] = 0;
    for (i = 0; i < V; ++i) voff[i+1] = voff[i] + vcount[i];
    {
        size_t *cursor = (size_t *)calloc(V, sizeof(size_t));
        if (!cursor) { rc = GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
        for (i = 0; i < E; ++i) {
            vneighbors[voff[edges[i].a] + cursor[edges[i].a]++] = edges[i].b;
            vneighbors[voff[edges[i].b] + cursor[edges[i].b]++] = edges[i].a;
        }
        free(cursor);
    }
    memcpy(rest, input->vertices, V*sizeof(geomlib_vec3));
    for (i = 0; i < V; ++i) {
        curx[i] = rest[i].x; cury[i] = rest[i].y; curz[i] = rest[i].z;
    }
    for (i = 0; i < handle_count; ++i) {
        if (handles[i].vertex >= V) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
        curx[handles[i].vertex] = handles[i].target.x;
        cury[handles[i].vertex] = handles[i].target.y;
        curz[handles[i].vertex] = handles[i].target.z;
    }

    /* Each edge contributes 3 components * 2 rows * 2 coefficients = 12
     * entries; each handle contributes 3 diagonal entries. */
    rows = (uint32_t *)calloc(12*E + 3*handle_count + 1, sizeof(uint32_t));
    cols = (uint32_t *)calloc(12*E + 3*handle_count + 1, sizeof(uint32_t));
    vals = (double *)calloc(12*E + 3*handle_count + 1, sizeof(double));
    solx = (double *)calloc(V, sizeof(double));
    soly = (double *)calloc(V, sizeof(double));
    solz = (double *)calloc(V, sizeof(double));
    if (!rows || !cols || !vals || !solx || !soly || !solz) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }
    sp = geomlib_sparse_solve_params_default();
    sp.max_iterations = 5000;
    sp.tolerance = 1e-6f;

    for (iter = 0; iter < (size_t)iterations; ++iter) {
        size_t count = 0, j;
        /* Local step: fit rotations from rest edges to current edges. */
        for (i = 0; i < V; ++i) {
            size_t n = voff[i+1]-voff[i];
            if (n == 0) { rot[i*9]=rot[i*9+4]=rot[i*9+8]=1.0; continue; }
            {
                geomlib_vec3 *aa = (geomlib_vec3 *)malloc(n*sizeof(geomlib_vec3));
                geomlib_vec3 *bb = (geomlib_vec3 *)malloc(n*sizeof(geomlib_vec3));
                if (!aa || !bb) { free(aa); free(bb); rc = GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
                for (j = 0; j < n; ++j) {
                    uint32_t nb = vneighbors[voff[i]+j];
                    aa[j] = geomlib_vec3_sub(rest[nb], rest[i]);
                    bb[j] = geomlib_vec3_create(
                        (float)(curx[nb]-curx[i]), (float)(cury[nb]-cury[i]),
                        (float)(curz[nb]-curz[i]));
                }
                fit_rotation(aa, bb, n, rot + i*9);
                free(aa); free(bb);
            }
        }
        /* Global step: solve for positions minimizing the ARAP residual.
         * Each edge contributes, per coordinate k:
         *   p_c - p_a = (R_a d)[k],   p_a - p_c = -(R_c d)[k],
         * with d = rest_c - rest_a. */
        {
            size_t nrows_eq = 0;
            double *bx = (double *)calloc(6*E + 3*handle_count + 1, sizeof(double));
            double *by = (double *)calloc(6*E + 3*handle_count + 1, sizeof(double));
            double *bz = (double *)calloc(6*E + 3*handle_count + 1, sizeof(double));
            if (!bx || !by || !bz) { free(bx); free(by); free(bz); rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
            for (i = 0; i < E; ++i) {
                uint32_t a = edges[i].a, c = edges[i].b;
                geomlib_vec3 dv = geomlib_vec3_sub(rest[c], rest[a]);
                double d3[3] = {dv.x, dv.y, dv.z};
                double e[2][3];
                int side, m;
                for (side = 0; side < 2; ++side) {
                    uint32_t center = side ? c : a;
                    const double *r = rot + center*9;
                    for (m = 0; m < 3; ++m) {
                        e[side][m] = r[m]*d3[0] + r[3+m]*d3[1] + r[6+m]*d3[2];
                    }
                }
                for (m = 0; m < 3; ++m) {
                    double *rhs = (m == 0) ? bx : (m == 1) ? by : bz;
                    /* p_c - p_a = e0 */
                    rows[count]=nrows_eq; cols[count]=a; vals[count++]=-1.0;
                    rows[count]=nrows_eq; cols[count]=c; vals[count++]= 1.0;
                    rhs[nrows_eq++] = e[0][m];
                    /* p_a - p_c = -e1 */
                    rows[count]=nrows_eq; cols[count]=a; vals[count++]= 1.0;
                    rows[count]=nrows_eq; cols[count]=c; vals[count++]=-1.0;
                    rhs[nrows_eq++] = -e[1][m];
                }
            }
            for (i = 0; i < handle_count; ++i) {
                uint32_t v = handles[i].vertex;
                int m;
                double target[3] = {handles[i].target.x, handles[i].target.y,
                                    handles[i].target.z};
                for (m = 0; m < 3; ++m) {
                    double *rhs = (m == 0) ? bx : (m == 1) ? by : bz;
                    rows[count]=nrows_eq; cols[count]=v; vals[count++]=1000.0;
                    rhs[nrows_eq++] = 1000.0 * target[m];
                }
            }
            {
                geomlib_spmat A = {0};
                rc = geomlib_spmat_from_triplets(rows, cols, vals, count,
                                                 nrows_eq, V, &A);
                if (rc != GEOMLIB_SPARSE_SUCCESS) { free(bx); free(by); free(bz); goto cleanup; }
                if (geomlib_solve_cgls(&A, bx, solx, &sp) != GEOMLIB_SPARSE_SUCCESS ||
                    geomlib_solve_cgls(&A, by, soly, &sp) != GEOMLIB_SPARSE_SUCCESS ||
                    geomlib_solve_cgls(&A, bz, solz, &sp) != GEOMLIB_SPARSE_SUCCESS) {
                    geomlib_spmat_free(&A); free(bx); free(by); free(bz);
                    rc = GEOMLIB_MESH_ERROR_INVALID;
                    goto cleanup;
                }
                geomlib_spmat_free(&A);
            }
            free(bx); free(by); free(bz);
        }
        for (i = 0; i < V; ++i) {
            curx[i] = solx[i]; cury[i] = soly[i]; curz[i] = solz[i];
        }
        for (i = 0; i < handle_count; ++i) {
            uint32_t v = handles[i].vertex;
            curx[v]=handles[i].target.x;
            cury[v]=handles[i].target.y;
            curz[v]=handles[i].target.z;
        }
    }

    {
        geomlib_trimesh result;
        rc = geomlib_trimesh_init(&result, V, T);
        if (rc != GEOMLIB_MESH_SUCCESS) goto cleanup;
        for (i = 0; i < V; ++i)
            geomlib_trimesh_add_vertex(&result,
                geomlib_vec3_create((float)curx[i], (float)cury[i], (float)curz[i]),
                NULL, NULL);
        for (i = 0; i < T; ++i)
            geomlib_trimesh_add_triangle(&result, input->indices[i*3],
                input->indices[i*3+1], input->indices[i*3+2]);
        *output = result;
    }

cleanup:
    free(edges); free(vcount); free(voff); free(vneighbors);
    free(rot); free(rest); free(curx); free(cury); free(curz);
    free(rows); free(cols); free(vals);
    free(solx); free(soly); free(solz);
    return rc;
}
