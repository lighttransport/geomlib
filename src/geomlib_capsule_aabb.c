/*
 * geomlib - Capsule / AABB overlap (see geomlib_capsule_aabb.h)
 *
 * Clean-room segment-box distance via endpoint projections and box edges.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_capsule_aabb.h"
#include "geomlib/geomlib_segment_distance.h"

#include <float.h>
#include <math.h>

static float clampv(float v, float lo, float hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static float point_box_dist(geomlib_vec3 q, geomlib_vec3 lo, geomlib_vec3 hi)
{
    geomlib_vec3 c = geomlib_vec3_create(
        clampv(q.x, lo.x, hi.x), clampv(q.y, lo.y, hi.y),
        clampv(q.z, lo.z, hi.z));
    return geomlib_vec3_length(geomlib_vec3_sub(q, c));
}

static int segment_hits_box(geomlib_vec3 a, geomlib_vec3 b,
                            geomlib_vec3 lo, geomlib_vec3 hi)
{
    float av[3] = {a.x, a.y, a.z}, bv[3] = {b.x, b.y, b.z};
    float l[3] = {lo.x, lo.y, lo.z}, h[3] = {hi.x, hi.y, hi.z};
    float tmin = 0.0f, tmax = 1.0f;
    int i;
    for (i = 0; i < 3; ++i) {
        float d = bv[i] - av[i], t0, t1;
        if (fabsf(d) < 1e-12f) {
            if (av[i] < l[i] || av[i] > h[i]) return 0;
            continue;
        }
        t0 = (l[i] - av[i]) / d;
        t1 = (h[i] - av[i]) / d;
        if (t0 > t1) { float t = t0; t0 = t1; t1 = t; }
        if (t0 > tmin) tmin = t0;
        if (t1 < tmax) tmax = t1;
        if (tmin > tmax) return 0;
    }
    return 1;
}

static geomlib_vec3 box_vertex(geomlib_vec3 lo, geomlib_vec3 hi, int code)
{
    return geomlib_vec3_create(
        (code & 1) ? hi.x : lo.x,
        (code & 2) ? hi.y : lo.y,
        (code & 4) ? hi.z : lo.z);
}

float geomlib_segment_aabb_distance(geomlib_vec3 a, geomlib_vec3 b,
                                    geomlib_vec3 lo, geomlib_vec3 hi)
{
    float best, da, db;
    int code, bit;

    if (segment_hits_box(a, b, lo, hi)) return 0.0f;
    da = point_box_dist(a, lo, hi);
    db = point_box_dist(b, lo, hi);
    best = da < db ? da : db;
    for (code = 0; code < 8; ++code)
        for (bit = 0; bit < 3; ++bit) {
            geomlib_vec3 p, q;
            geomlib_segment_distance sd;
            if (code & (1 << bit)) continue;
            p = box_vertex(lo, hi, code);
            q = box_vertex(lo, hi, code | (1 << bit));
            if (geomlib_segment_closest_points(a, b, p, q, &sd) ==
                GEOMLIB_MESH_SUCCESS && sd.distance < best)
                best = sd.distance;
        }
    return best;
}

int geomlib_capsule_aabb_overlap(geomlib_vec3 a, geomlib_vec3 b, float r,
                                 geomlib_vec3 lo, geomlib_vec3 hi,
                                 geomlib_capsule_aabb_result *out)
{
    geomlib_capsule_aabb_result res;
    float d;

    res.intersects = 0;
    res.distance = 0.0f;
    if (!out || !(r >= 0.0f)) return 0;
    d = geomlib_segment_aabb_distance(a, b, lo, hi);
    res.distance = d - r;
    res.intersects = (res.distance <= 0.0f);
    *out = res;
    return res.intersects;
}
