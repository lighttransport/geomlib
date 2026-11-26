/*
 * geomlib - Cotangent Laplacian (see geomlib_laplacian.h)
 *
 * Clean-room cotangent stiffness + lumped mass assembly.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_laplacian.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_laplacian geomlib_laplacian_create(void)
{
    geomlib_laplacian l;
    memset(&l,0,sizeof(l));
    return l;
}

void geomlib_laplacian_free(geomlib_laplacian *l)
{
    if (!l) return;
    geomlib_spmat_free(&l->stiffness);
    free(l->mass);
    memset(l,0,sizeof(*l));
}

int geomlib_build_laplacian(const geomlib_trimesh *mesh, geomlib_laplacian *out)
{
    uint32_t *rows=NULL,*cols=NULL;
    double *vals=NULL,*mass=NULL,*diag=NULL;
    size_t n,ne=0,i;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out||!mesh||!mesh->vertices||!mesh->indices||!mesh->triangle_count)
        return rc;
    n=mesh->vertex_count;
    rows=(uint32_t*)malloc((6*mesh->triangle_count+n)*sizeof(uint32_t));
    cols=(uint32_t*)malloc((6*mesh->triangle_count+n)*sizeof(uint32_t));
    vals=(double*)malloc((6*mesh->triangle_count+n)*sizeof(double));
    mass=(double*)calloc(n,sizeof(double));
    diag=(double*)calloc(n,sizeof(double));
    if (!rows||!cols||!vals||!mass||!diag) goto cleanup;
    memset(out,0,sizeof(*out));

    for (i=0;i<mesh->triangle_count;++i) {
        uint32_t v[3]={mesh->indices[i*3],mesh->indices[i*3+1],mesh->indices[i*3+2]};
        geomlib_vec3d p[3];
        geomlib_vec3d normal,cr;
        double area;
        int k;
        if (v[0]>=n||v[1]>=n||v[2]>=n) goto cleanup;
        for (k=0;k<3;++k) {
            p[k].x=mesh->vertices[v[k]].x;
            p[k].y=mesh->vertices[v[k]].y;
            p[k].z=mesh->vertices[v[k]].z;
        }
        normal=geomlib_vec3d_cross(geomlib_vec3d_sub(p[1],p[0]),
                                   geomlib_vec3d_sub(p[2],p[0]));
        area=0.5*geomlib_vec3d_length(normal);
        for (k=0;k<3;++k) mass[v[k]]+=area/3.0;
        for (k=0;k<3;++k) {
            int a=k,b=(k+1)%3,c=(k+2)%3;
            geomlib_vec3d u=geomlib_vec3d_sub(p[c],p[a]);
            geomlib_vec3d v2=geomlib_vec3d_sub(p[c],p[b]);
            cr=geomlib_vec3d_cross(u,v2);
            double l=geomlib_vec3d_length(cr);
            if (!(l>1e-14)) continue;
            double w=0.5*geomlib_vec3d_dot(u,v2)/l;
            rows[ne]=v[a]; cols[ne]=v[b]; vals[ne]=w; diag[v[a]]-=w; ne++;
            rows[ne]=v[b]; cols[ne]=v[a]; vals[ne]=w; diag[v[b]]-=w; ne++;
        }
    }
    for (i=0;i<n;++i) {
        rows[ne]=(uint32_t)i; cols[ne]=(uint32_t)i; vals[ne]=diag[i]; ne++;
    }
    rc=geomlib_spmat_from_triplets(rows,cols,vals,ne,n,n,&out->stiffness);
    if (rc!=GEOMLIB_SPARSE_SUCCESS) goto cleanup;
    out->mass=mass; mass=NULL;
    out->vertex_count=n;
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    free(rows); free(cols); free(vals); free(mass); free(diag);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_laplacian_free(out);
    return rc;
}

int geomlib_laplacian_apply(const geomlib_laplacian *l, const double *x,
                            double scale, double *y)
{
    size_t i,k;
    if (!l||!x||!y||!l->stiffness.values) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<l->vertex_count;++i) {
        double sum=0.0;
        for (k=l->stiffness.rowptr[i];k<l->stiffness.rowptr[i+1];++k)
            sum+=l->stiffness.values[k]*x[l->stiffness.colidx[k]];
        y[i]=sum+scale*l->mass[i]*x[i];
    }
    return GEOMLIB_MESH_SUCCESS;
}
