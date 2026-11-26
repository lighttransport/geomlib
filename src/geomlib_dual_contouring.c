/*
 * geomlib - Dual contouring (see geomlib_dual_contouring.h)
 *
 * Clean-room implementation of dual contouring (Ju, Losasso, Schaefer,
 * Hanrahan, "Dual Contouring of Hermite Data", SIGGRAPH 2002) on a uniform
 * grid: hermite data (position + normal) comes from the SDF zero crossings;
 * each crossing cell gets one vertex minimizing the sum of plane quadrics;
 * quads span the four cells around each crossing grid edge.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_dual_contouring.h"
#include "geomlib/geomlib_repair.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double a[3][3];
    double b[3];
    double w;
} quadric;

static void quadric_add(quadric *q, double nx, double ny, double nz,
                        double px, double py, double pz, double weight)
{
    double d = nx*px + ny*py + nz*pz;
    double v[3] = {nx, ny, nz};
    int i, j;
    for (i = 0; i < 3; ++i) {
        q->b[i] += weight * d * v[i];
        for (j = 0; j < 3; ++j) q->a[i][j] += weight * v[i]*v[j];
    }
    q->w += weight;
}

/* Solve the 3x3 normal equations with partial pivoting. */
static geomlib_vec3 quadric_solve(const quadric *q, geomlib_vec3 fallback)
{
    double m[3][4];
    int i, j, k;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) m[i][j] = q->a[i][j] + (i==j ? 1e-8*q->w : 0.0);
        m[i][3] = q->b[i];
    }
    for (k = 0; k < 3; ++k) {
        int pivot = k;
        for (i = k+1; i < 3; ++i) if (fabs(m[i][k]) > fabs(m[pivot][k])) pivot = i;
        if (fabs(m[pivot][k]) < 1e-12) return fallback;
        if (pivot != k) {
            for (j = 0; j < 4; ++j) { double t=m[k][j]; m[k][j]=m[pivot][j]; m[pivot][j]=t; }
        }
        for (i = k+1; i < 3; ++i) {
            double f = m[i][k]/m[k][k];
            for (j = k; j < 4; ++j) m[i][j] -= f*m[k][j];
        }
    }
    {
        double x[3];
        for (k = 2; k >= 0; --k) {
            double s = m[k][3];
            for (j = k+1; j < 3; ++j) s -= m[k][j]*x[j];
            x[k] = s/m[k][k];
        }
        return geomlib_vec3_create((float)x[0], (float)x[1], (float)x[2]);
    }
}

typedef struct {
    geomlib_vec3 *v;
    size_t count, cap;
    uint32_t *idx;
    size_t icount, icap;
} soup;

static void add_tri(soup *s, uint32_t a, uint32_t b, uint32_t c, geomlib_vec3 outward)
{
    geomlib_vec3 n = geomlib_vec3_cross(
        geomlib_vec3_sub(s->v[b], s->v[a]), geomlib_vec3_sub(s->v[c], s->v[a]));
    if (geomlib_vec3_dot(n, outward) < 0.0f) {
        uint32_t t = b; b = c; c = t;
    }
    if (s->icount + 3 > s->icap) {
        size_t nc = s->icap ? s->icap*2 : 4096;
        uint32_t *ni = (uint32_t *)realloc(s->idx, nc*sizeof(uint32_t));
        if (!ni) return;
        s->idx = ni;
        s->icap = nc;
    }
    s->idx[s->icount++] = a;
    s->idx[s->icount++] = b;
    s->idx[s->icount++] = c;
}

static double sample(const geomlib_sdf *s, size_t x, size_t y, size_t z)
{
    return s->values[x + s->resolution[0]*(y + s->resolution[1]*z)];
}

static geomlib_vec3 gradient(const geomlib_sdf *s, size_t x, size_t y, size_t z)
{
    size_t nx = s->resolution[0], ny = s->resolution[1], nz = s->resolution[2];
    geomlib_vec3 g;
    g.x = (float)(sample(s, (x+1<nx)?x+1:nx-1, y, z) - sample(s, x>0?x-1:0, y, z));
    g.y = (float)(sample(s, x, (y+1<ny)?y+1:ny-1, z) - sample(s, x, y>0?y-1:0, z));
    g.z = (float)(sample(s, x, y, (z+1<nz)?z+1:nz-1) - sample(s, x, y, z>0?z-1:0));
    return geomlib_vec3_normalize(g);
}

int geomlib_sdf_dual_contouring(const geomlib_sdf *sdf, geomlib_trimesh *out)
{
    size_t nx, ny, nz, x, y, z, i;
    int32_t *cell_vertex;
    geomlib_vec3 *cell_pos;
    soup s;

    if (!sdf || !out || !sdf->values) return GEOMLIB_MESH_ERROR_INVALID;
    nx = sdf->resolution[0]; ny = sdf->resolution[1]; nz = sdf->resolution[2];
    if (nx < 3 || ny < 3 || nz < 3) return GEOMLIB_MESH_ERROR_INVALID;
    {
        size_t cells = (nx-1)*(ny-1)*(nz-1);
        cell_vertex = (int32_t *)malloc(cells*sizeof(int32_t));
        cell_pos = (geomlib_vec3 *)malloc(cells*sizeof(geomlib_vec3));
        if (!cell_vertex || !cell_pos) {
            free(cell_vertex); free(cell_pos);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    for (i = 0; i < (nx-1)*(ny-1)*(nz-1); ++i) cell_vertex[i] = -1;
    memset(&s, 0, sizeof(s));
#define CELL_ID(x,y,z) ((x) + (nx-1)*((y) + (ny-1)*(z)))

    /* Place one vertex per crossing cell. */
    for (z = 0; z+1 < nz; ++z) {
        for (y = 0; y+1 < ny; ++y) {
            for (x = 0; x+1 < nx; ++x) {
                quadric q;
                double sum[3] = {0,0,0};
                size_t hits = 0;
                memset(&q, 0, sizeof(q));
                /* Check the 12 cell edges for sign changes. */
                for (i = 0; i < 12; ++i) {
                    size_t ax, ay, az, bx, by, bz;
                    static const int offs[12][2][3] = {
                        {{0,0,0},{1,0,0}}, {{0,1,0},{1,1,0}}, {{0,0,1},{1,0,1}}, {{0,1,1},{1,1,1}},
                        {{0,0,0},{0,1,0}}, {{1,0,0},{1,1,0}}, {{0,0,1},{0,1,1}}, {{1,0,1},{1,1,1}},
                        {{0,0,0},{0,0,1}}, {{1,0,0},{1,0,1}}, {{0,1,0},{0,1,1}}, {{1,1,0},{1,1,1}}
                    };
                    double va, vb, t;
                    geomlib_vec3 pa, pb, hit, n;
                    ax=x+offs[i][0][0]; ay=y+offs[i][0][1]; az=z+offs[i][0][2];
                    bx=x+offs[i][1][0]; by=y+offs[i][1][1]; bz=z+offs[i][1][2];
                    va = sample(sdf, ax,ay,az);
                    vb = sample(sdf, bx,by,bz);
                    if ((va < 0) == (vb < 0)) continue;
                    t = va / (va - vb);
                    pa = geomlib_vec3_create(sdf->origin.x+(float)ax*sdf->voxel_size,
                                             sdf->origin.y+(float)ay*sdf->voxel_size,
                                             sdf->origin.z+(float)az*sdf->voxel_size);
                    pb = geomlib_vec3_create(sdf->origin.x+(float)bx*sdf->voxel_size,
                                             sdf->origin.y+(float)by*sdf->voxel_size,
                                             sdf->origin.z+(float)bz*sdf->voxel_size);
                    hit = geomlib_vec3_add(geomlib_vec3_mul(pa,(float)(1-t)),
                                           geomlib_vec3_mul(pb,(float)t));
                    n = gradient(sdf, ax,ay,az);
                    quadric_add(&q, n.x, n.y, n.z, hit.x, hit.y, hit.z, 1.0);
                    sum[0] += hit.x; sum[1] += hit.y; sum[2] += hit.z;
                    hits++;
                }
                if (hits == 0) continue;
                {
                    geomlib_vec3 fallback = geomlib_vec3_create(
                        (float)(sum[0]/hits), (float)(sum[1]/hits), (float)(sum[2]/hits));
                    geomlib_vec3 p = quadric_solve(&q, fallback);
                    size_t cid = CELL_ID(x,y,z);
                    float md = sdf->voxel_size;
                    if (p.x < sdf->origin.x-md || p.x > sdf->origin.x+nx*sdf->voxel_size+md ||
                        p.y < sdf->origin.y-md || p.y > sdf->origin.y+ny*sdf->voxel_size+md ||
                        p.z < sdf->origin.z-md || p.z > sdf->origin.z+nz*sdf->voxel_size+md)
                        p = fallback;
                    if (s.count == s.cap) {
                        size_t nc = s.cap ? s.cap*2 : 4096;
                        geomlib_vec3 *nv = (geomlib_vec3 *)realloc(s.v, nc*sizeof(geomlib_vec3));
                        if (!nv) { free(s.v); free(s.idx); free(cell_vertex); free(cell_pos); return GEOMLIB_MESH_ERROR_MEMORY; }
                        s.v = nv; s.cap = nc;
                    }
                    s.v[s.count] = p;
                    cell_vertex[cid] = (int32_t)s.count++;
                    cell_pos[cid] = p;
                }
            }
        }
    }

    /* Emit quads around each crossing grid edge (x, y, z directed). */
    for (z = 1; z+1 < nz; ++z) {
        for (y = 1; y+1 < ny; ++y) {
            for (x = 1; x+1 < nx; ++x) {
                size_t axis;
                for (axis = 0; axis < 3; ++axis) {
                    size_t gx=x, gy=y, gz=z;
                    uint32_t ids[4];
                    size_t n = 0, a, b;
                    /* Cells around a primal edge: the cell index along the
                     * edge direction is fixed; the other two vary. */
                    static const int dirs[3][4][3] = {
                        {{0,-1,-1},{0,-1,0},{0,0,0},{0,0,-1}},
                        {{-1,0,-1},{0,0,-1},{0,0,0},{-1,0,0}},
                        {{-1,-1,0},{0,-1,0},{0,0,0},{-1,0,0}}
                    };
                    /* Cells around the grid edge at (x,y,z) along `axis`. */
                    for (a = 0; a < 4; ++a) {
                        int dx = dirs[axis][a][0];
                        int dy = dirs[axis][a][1];
                        int dz = dirs[axis][a][2];
                        size_t cx = (size_t)((int)gx + dx);
                        size_t cy = (size_t)((int)gy + dy);
                        size_t cz = (size_t)((int)gz + dz);
                        int32_t cv;
                        if (cx == 0 || cy == 0 || cz == 0 ||
                            cx >= nx-1 || cy >= ny-1 || cz >= nz-1) continue;
                        cv = cell_vertex[CELL_ID(cx,cy,cz)];
                        if (cv < 0) continue;
                        ids[n++] = (uint32_t)cv;
                    }
                    if (n != 4) continue;
                    /* Skip duplicate vertices (degenerate quad). */
                    for (a = 0; a < 4; ++a)
                        for (b = a+1; b < 4; ++b)
                            if (ids[a] == ids[b]) n = 0;
                    if (n != 4) continue;
                    {
                        geomlib_vec3 grid_point = {
                            sdf->origin.x + (float)x*sdf->voxel_size,
                            sdf->origin.y + (float)y*sdf->voxel_size,
                            sdf->origin.z + (float)z*sdf->voxel_size
                        };
                        geomlib_vec3 avg = {0,0,0};
                        geomlib_vec3 outward = gradient(sdf, (uint32_t)x,
                                                        (uint32_t)y, (uint32_t)z);
                        (void)avg;
                        (void)grid_point;
                        add_tri(&s, ids[0], ids[1], ids[2], outward);
                        add_tri(&s, ids[0], ids[2], ids[3], outward);
                    }
                }
            }
        }
    }
#undef CELL_ID

    if (s.icount == 0) {
        geomlib_trimesh_init(out, 0, 0);
        free(s.v); free(s.idx); free(cell_vertex); free(cell_pos);
        return GEOMLIB_MESH_SUCCESS;
    }
    {
        int rc = geomlib_soup_to_mesh(s.v, s.count, s.idx, s.icount,
                                      sdf->voxel_size*1e-3f, out);
        free(s.v); free(s.idx); free(cell_vertex); free(cell_pos);
        return rc;
    }
}
