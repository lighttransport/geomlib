/*
 * geomlib - Solid angle (see geomlib_solid_angle.h)
 *
 * Clean-room Van Oosterom-Strackee formula. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_solid_angle.h"

#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void to_double(geomlib_vec3 v, double *out)
{
    out[0] = (double)v.x;
    out[1] = (double)v.y;
    out[2] = (double)v.z;
}

double geomlib_solid_angle_triangle(geomlib_vec3 p, geomlib_vec3 a,
                                    geomlib_vec3 b, geomlib_vec3 c)
{
    double va[3], vb[3], vc[3], na, nb, nc, num, den;
    to_double(geomlib_vec3_sub(a, p), va);
    to_double(geomlib_vec3_sub(b, p), vb);
    to_double(geomlib_vec3_sub(c, p), vc);
    na = sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
    nb = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);
    nc = sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);
    if (!(na > 1e-18 && nb > 1e-18 && nc > 1e-18)) return 0.0;

    /* num = a . (b x c); den = |a||b||c| + (a.b)|c| + (b.c)|a| + (c.a)|b|. */
    {
        double cx = vb[1] * vc[2] - vb[2] * vc[1];
        double cy = vb[2] * vc[0] - vb[0] * vc[2];
        double cz = vb[0] * vc[1] - vb[1] * vc[0];
        double ab = va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
        double bc = vb[0] * vc[0] + vb[1] * vc[1] + vb[2] * vc[2];
        double ca = vc[0] * va[0] + vc[1] * va[1] + vc[2] * va[2];
        num = va[0] * cx + va[1] * cy + va[2] * cz;
        den = na * nb * nc + ab * nc + bc * na + ca * nb;
    }
    return 2.0 * atan2(num, den);
}

double geomlib_solid_angle_polygon(const geomlib_vec3 *verts, size_t n,
                                   geomlib_vec3 p)
{
    size_t i;
    double total = 0.0;
    if (!verts || n < 3) return 0.0;
    for (i = 1; i + 1 < n; ++i)
        total += geomlib_solid_angle_triangle(p, verts[0], verts[i], verts[i + 1]);
    return total;
}

double geomlib_solid_angle_mesh(const geomlib_trimesh *mesh, geomlib_vec3 p)
{
    size_t i;
    double total = 0.0;
    if (!mesh || !mesh->indices || !mesh->vertices) return 0.0;
    for (i = 0; i < mesh->triangle_count; ++i) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];
        if (i0 >= mesh->vertex_count || i1 >= mesh->vertex_count ||
            i2 >= mesh->vertex_count) continue;
        total += geomlib_solid_angle_triangle(
            p, mesh->vertices[i0], mesh->vertices[i1], mesh->vertices[i2]);
    }
    return total;
}

double geomlib_winding_number_mesh(const geomlib_trimesh *mesh, geomlib_vec3 p)
{
    return geomlib_solid_angle_mesh(mesh, p) / (4.0 * M_PI);
}
