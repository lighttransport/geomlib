/*
 * geomlib - Triangle-mesh voxelization (see geomlib_voxelize.h)
 *
 * Clean-room implementation. Surface voxels are found with the 13-axis
 * separating-axis test between a triangle and an AABB (Akenine-Moller).
 * Interior voxels use ray parity against a BVH along +X.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_voxelize.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void geomlib_voxel_grid_free(geomlib_voxel_grid *g)
{
    if (!g) return;
    free(g->surface);
    free(g->volume);
    memset(g, 0, sizeof(*g));
}

size_t geomlib_voxel_grid_index(const geomlib_voxel_grid *g,
                                uint32_t x, uint32_t y, uint32_t z)
{
    if (!g || !g->surface) return (size_t)-1;
    if (x >= g->nx || y >= g->ny || z >= g->nz) return (size_t)-1;
    return (size_t)x + (size_t)g->nx*((size_t)y + (size_t)g->ny*z);
}

typedef struct { double x, y, z; } vec3d;

static vec3d vsub(vec3d a, vec3d b) { vec3d r={a.x-b.x,a.y-b.y,a.z-b.z}; return r; }
static vec3d cross(vec3d a, vec3d b)
{
    vec3d r;
    r.x = a.y*b.z - a.z*b.y;
    r.y = a.z*b.x - a.x*b.z;
    r.z = a.x*b.y - a.y*b.x;
    return r;
}
static double dot(vec3d a, vec3d b) { return a.x*b.x+a.y*b.y+a.z*b.z; }

/* Does the axis L separate the triangle from the box? */
static int separated(vec3d L, const vec3d t[3], vec3d c,
                     double hx, double hy, double hz)
{
    double p0=dot(L,t[0]), p1=dot(L,t[1]), p2=dot(L,t[2]);
    double tmin = p0<p1 ? (p0<p2?p0:p2) : (p1<p2?p1:p2);
    double tmax = p0>p1 ? (p0>p2?p0:p2) : (p1>p2?p1:p2);
    double cr = hx*fabs(L.x) + hy*fabs(L.y) + hz*fabs(L.z);
    double cc = dot(L,c);
    return tmax < cc-cr || tmin > cc+cr;
}

/* 13-axis SAT triangle/AABB overlap. */
static int tri_box_overlap(const geomlib_vec3 *vp, vec3d c,
                           double hx, double hy, double hz)
{
    vec3d t[3], e0, e1, e2, n, axis;
    int i;
    for (i = 0; i < 3; ++i) {
        t[i].x = vp[i].x; t[i].y = vp[i].y; t[i].z = vp[i].z;
    }
    /* Box face normals. */
    axis.x=1; axis.y=0; axis.z=0; if (separated(axis,t,c,hx,hy,hz)) return 0;
    axis.x=0; axis.y=1; axis.z=0; if (separated(axis,t,c,hx,hy,hz)) return 0;
    axis.x=0; axis.y=0; axis.z=1; if (separated(axis,t,c,hx,hy,hz)) return 0;
    e0 = vsub(t[1],t[0]); e1 = vsub(t[2],t[1]); e2 = vsub(t[0],t[2]);
    /* Triangle normal. */
    n = cross(e0,e1);
    if (dot(n,n) > 1e-20 && separated(n,t,c,hx,hy,hz)) return 0;
    /* Nine cross-product axes. */
    {
        const double u[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        const vec3d edges[3] = {e0,e1,e2};
        int a, b;
        for (a = 0; a < 3; ++a)
            for (b = 0; b < 3; ++b) {
                vec3d U = {u[a][0],u[a][1],u[a][2]};
                axis = cross(U,edges[b]);
                if (dot(axis,axis) < 1e-20) continue;
                if (separated(axis,t,c,hx,hy,hz)) return 0;
            }
    }
    return 1;
}

geomlib_voxel_grid geomlib_voxelize_mesh(const geomlib_trimesh *mesh,
                                         uint32_t nx, uint32_t ny, uint32_t nz)
{
    geomlib_voxel_grid g;
    geomlib_aabb box;
    geomlib_vec3 size, ext;
    double hx, hy, hz;
    size_t total, i;
    uint32_t x, y, z;
    geomlib_bvh bvh;
    int bvh_ok = 0;

    memset(&g, 0, sizeof(g));
    if (!mesh || !mesh->vertices || !mesh->indices ||
        !nx || !ny || !nz || mesh->triangle_count == 0) return g;
    total = (size_t)nx*ny*nz;
    if (total > 64u*1024u*1024u) return g;

    g.surface = (uint8_t *)calloc(total, 1);
    g.volume = (uint8_t *)calloc(total, 1);
    if (!g.surface || !g.volume) { geomlib_voxel_grid_free(&g); return g; }

    box = geomlib_trimesh_compute_aabb(mesh);
    size = geomlib_aabb_size(box);
    ext = geomlib_vec3_create(fmaxf(size.x, 1e-4f),
                              fmaxf(size.y, 1e-4f),
                              fmaxf(size.z, 1e-4f));
    g.voxel_size = 1.10f*fmaxf(ext.x, fmaxf(ext.y, ext.z)) /
                   (float)(nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz));
    g.origin = geomlib_vec3_sub(geomlib_aabb_center(box),
        geomlib_vec3_mul(ext, 0.5f + 0.05f));
    g.nx = nx; g.ny = ny; g.nz = nz;
    hx = hy = hz = 0.5*(double)g.voxel_size;

    /* Surface: conservative SAT overlap per triangle over its voxel range. */
    for (i = 0; i < mesh->triangle_count; ++i) {
        const geomlib_vec3 *tri;
        uint32_t lo[3], hi[3];
        int axis;
        geomlib_vec3 vs[3];
        for (int k = 0; k < 3; ++k) {
            uint32_t id = mesh->indices[i*3+k];
            if (id >= mesh->vertex_count) { vs[k] = geomlib_vec3_create(0,0,0); }
            else vs[k] = mesh->vertices[id];
        }
        tri = vs;
        for (axis = 0; axis < 3; ++axis) {
            float f0 = (axis==0)?vs[0].x:((axis==1)?vs[0].y:vs[0].z);
            float f1 = (axis==0)?vs[1].x:((axis==1)?vs[1].y:vs[1].z);
            float f2 = (axis==0)?vs[2].x:((axis==1)?vs[2].y:vs[2].z);
            float o = (axis==0)?g.origin.x:((axis==1)?g.origin.y:g.origin.z);
            float mn = f0<f1?(f0<f2?f0:f2):(f1<f2?f1:f2);
            float mx = f0>f1?(f0>f2?f0:f2):(f1>f2?f1:f2);
            int l = (int)floorf((mn-o)/g.voxel_size)-1;
            int h = (int)ceilf((mx-o)/g.voxel_size)+1;
            uint32_t n = (axis==0)?nx:((axis==1)?ny:nz);
            lo[axis] = l < 0 ? 0 : (uint32_t)l;
            hi[axis] = h >= (int)n ? n-1 : (uint32_t)h;
        }
        for (z = lo[2]; z <= hi[2]; ++z)
            for (y = lo[1]; y <= hi[1]; ++y)
                for (x = lo[0]; x <= hi[0]; ++x) {
                    size_t id = geomlib_voxel_grid_index(&g, x, y, z);
                    geomlib_vec3 c;
                    vec3d ctr;
                    if (g.surface[id]) continue;
                    c = geomlib_vec3_create(
                        g.origin.x + ((float)x+0.5f)*g.voxel_size,
                        g.origin.y + ((float)y+0.5f)*g.voxel_size,
                        g.origin.z + ((float)z+0.5f)*g.voxel_size);
                    ctr.x=c.x; ctr.y=c.y; ctr.z=c.z;
                    if (tri_box_overlap(tri, ctr, hx, hy, hz)) g.surface[id] = 1;
                }
    }

    /* Volume: ray parity along +X. */
    if (geomlib_bvh_build(&bvh, mesh, GEOMLIB_BVH_QUALITY_MEDIUM) == GEOMLIB_MESH_SUCCESS)
        bvh_ok = 1;
    if (bvh_ok) {
        geomlib_vec3 dir = geomlib_vec3_create(1.0f, 0.0f, 0.0f);
        for (z = 0; z < nz; ++z)
            for (y = 0; y < ny; ++y)
                for (x = 0; x < nx; ++x) {
                    size_t id = geomlib_voxel_grid_index(&g, x, y, z);
                    geomlib_vec3 o = geomlib_vec3_create(
                        g.origin.x + ((float)x+0.5f)*g.voxel_size,
                        g.origin.y + ((float)y+0.5f)*g.voxel_size,
                        g.origin.z + ((float)z+0.5f)*g.voxel_size);
                    int crossings = 0, iter;
                    for (iter = 0; iter < 64; ++iter) {
                        geomlib_ray ray = geomlib_ray_create(o, dir);
                        geomlib_ray_hit hit;
                        ray.t_min = 1e-6f;
                        if (!geomlib_bvh_trace_ray(&bvh, &ray, &hit)) break;
                        crossings++;
                        o = geomlib_vec3_add(hit.position,
                            geomlib_vec3_create(1e-4f, 0.0f, 0.0f));
                    }
                    if (crossings & 1) g.volume[id] = 1;
                }
        geomlib_bvh_free(&bvh);
    }

    for (i = 0; i < total; ++i) {
        if (g.surface[i]) g.surface_count++;
        if (g.volume[i]) g.volume_count++;
    }
    return g;
}
