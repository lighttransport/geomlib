/*
 * geomlib - 2D polygon distance (see geomlib_polygon_distance2d.h)
 *
 * Clean-room closest-feature enumeration using existing point-in-polygon,
 * segment intersection and segment closest-point helpers. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_distance2d.h"
#include "geomlib/geomlib_point_in_polygon.h"
#include "geomlib/geomlib_segment2d.h"
#include "geomlib/geomlib_segment_distance.h"

#include <float.h>

static geomlib_vec3 to3(geomlib_vec2 p)
{
    return geomlib_vec3_create(p.x, p.y, 0.0f);
}

static geomlib_vec2 to2(geomlib_vec3 p)
{
    return geomlib_vec2_create(p.x, p.y);
}

static void consider(geomlib_polygon_distance2d_result *best,
                     geomlib_vec2 pa, geomlib_vec2 pb)
{
    geomlib_vec3 a = to3(pa), b = to3(pb);
    float d = geomlib_vec3_length(geomlib_vec3_sub(a, b));
    if (d < best->distance) {
        best->distance = d;
        best->point_a = pa;
        best->point_b = pb;
    }
}

int geomlib_polygon_distance2d(const geomlib_vec2 *a, size_t na,
                               const geomlib_vec2 *b, size_t nb,
                               geomlib_polygon_distance2d_result *out)
{
    geomlib_polygon_distance2d_result r;
    size_t i, j;

    r.distance = FLT_MAX;
    r.point_a = geomlib_vec2_create(0, 0);
    r.point_b = geomlib_vec2_create(0, 0);
    if (!out || !a || !b || na < 3 || nb < 3) return 0;

    /* Overlap: a vertex inside the other polygon, or crossing edges. */
    for (i = 0; i < na; ++i)
        if (geomlib_point_in_polygon(b, nb, a[i])) {
            r.distance = 0.0f;
            r.point_a = r.point_b = a[i];
            *out = r;
            return 1;
        }
    for (i = 0; i < nb; ++i)
        if (geomlib_point_in_polygon(a, na, b[i])) {
            r.distance = 0.0f;
            r.point_a = r.point_b = b[i];
            *out = r;
            return 1;
        }
    for (i = 0; i < na; ++i)
        for (j = 0; j < nb; ++j) {
            geomlib_seg2d_result s = geomlib_segment2d_intersect(
                a[i], a[(i + 1) % na], b[j], b[(j + 1) % nb]);
            if (s.kind != GEOMLIB_SEG2D_NONE) {
                r.distance = 0.0f;
                r.point_a = r.point_b = s.point;
                *out = r;
                return 1;
            }
        }

    /* For disjoint simple polygons the closest boundary pair is a vertex of
     * one polygon projected onto an edge of the other. */
    for (i = 0; i < na; ++i)
        for (j = 0; j < nb; ++j) {
            geomlib_segment_distance sd;
            if (geomlib_segment_closest_points(
                    to3(a[i]), to3(a[i]),
                    to3(b[j]), to3(b[(j + 1) % nb]), &sd) == GEOMLIB_MESH_SUCCESS)
                consider(&r, a[i], to2(sd.point_b));
        }
    for (i = 0; i < nb; ++i)
        for (j = 0; j < na; ++j) {
            geomlib_segment_distance sd;
            if (geomlib_segment_closest_points(
                    to3(b[i]), to3(b[i]),
                    to3(a[j]), to3(a[(j + 1) % na]), &sd) == GEOMLIB_MESH_SUCCESS)
                consider(&r, to2(sd.point_b), b[i]);
        }
    *out = r;
    return 1;
}
