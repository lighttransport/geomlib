/*
 * Wave 1 parameterization tests: cotangent harmonic, free-boundary LSCM and
 * ARAP UV on a triangulated disk, plus honest distortion statistics.
 */
#include <geomlib/geomlib_parameterization.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static geomlib_trimesh make_disk(size_t segments)
{
    geomlib_trimesh m;
    geomlib_trimesh_init(&m, segments + 1, segments);
    geomlib_trimesh_add_vertex(&m, geomlib_vec3_create(0, 0, 0), 0, 0);
    for (size_t i = 0; i < segments; ++i) {
        float a = 2.0f * 3.14159265f * (float)i / (float)segments;
        geomlib_trimesh_add_vertex(&m,
            geomlib_vec3_create(cosf(a), sinf(a), 0.02f * sinf(3*a)), 0, 0);
    }
    for (size_t i = 0; i < segments; ++i)
        geomlib_trimesh_add_triangle(&m, 0, (uint32_t)(i+1), (uint32_t)((i+1)%segments+1));
    return m;
}

static void check_uvs_valid(const geomlib_trimesh *m, const geomlib_uv_coords *uv,
                            const char *name)
{
    for (size_t i = 0; i < uv->vertex_count; ++i) {
        if (!std::isfinite(uv->uvs[i].x) || !std::isfinite(uv->uvs[i].y)) {
            std::printf("FAIL %s: non-finite UV at %zu\n", name, i);
            failures++;
            return;
        }
    }
    geomlib_param_stats stats;
    geomlib_compute_param_stats(m, uv, &stats);
    check(!stats.has_flips, name);
    check(stats.area_distortion >= 0.0f && std::isfinite(stats.area_distortion), "area distortion finite");
    check(stats.stretch_L2 > 0.0f && std::isfinite(stats.stretch_L2), "stretch L2 positive");
    check(stats.angle_distortion >= 0.0f, "angle distortion non-negative");
}

int main(void)
{
    geomlib_trimesh disk = make_disk(12);

    /* Cotangent harmonic with a circular boundary. */
    {
        geomlib_uv_coords uv = geomlib_uv_coords_create(disk.vertex_count, disk.triangle_count);
        check(geomlib_parameterize_harmonic(&disk, GEOMLIB_BOUNDARY_CIRCLE, &uv) ==
              GEOMLIB_MESH_SUCCESS, "harmonic rc");
        check_uvs_valid(&disk, &uv, "harmonic has no flips");
        check(std::fabs(uv.uvs[1].x - 1.0f) < 1e-5f && std::fabs(uv.uvs[1].y - 0.5f) < 1e-5f,
              "harmonic pins first boundary vertex to circle");
        geomlib_compute_param_distortion(&disk, &uv);
        geomlib_param_stats stats;
        geomlib_compute_param_stats(&disk, &uv, &stats);
        check(stats.area_distortion > 0.0f, "area distortion computed");
        geomlib_uv_coords_free(&uv);
    }

    /* Free-boundary LSCM: only two vertices pinned. */
    {
        geomlib_uv_coords uv = geomlib_uv_coords_create(disk.vertex_count, disk.triangle_count);
        check(geomlib_parameterize_lscm_free(&disk, 1, 7, &uv) == GEOMLIB_MESH_SUCCESS,
              "lscm free rc");
        check_uvs_valid(&disk, &uv, "lscm free has no flips");
        check(std::fabs(uv.uvs[1].x) < 1e-4f && std::fabs(uv.uvs[1].y) < 1e-4f,
              "lscm pins first vertex at origin");
        check(std::fabs(uv.uvs[7].x - 1.0f) < 1e-4f, "lscm pins second vertex at (1,0)");
        geomlib_uv_coords_free(&uv);
    }

    /* ARAP with boundary pinned to a circle. */
    {
        geomlib_uv_coords uv = geomlib_uv_coords_create(disk.vertex_count, disk.triangle_count);
        check(geomlib_parameterize_arap(&disk, GEOMLIB_BOUNDARY_CIRCLE, 4, &uv) ==
              GEOMLIB_MESH_SUCCESS, "arap rc");
        check_uvs_valid(&disk, &uv, "arap has no flips");
        check(std::fabs(uv.uvs[1].x - 1.0f) < 1e-3f, "arap keeps boundary pin");
        geomlib_uv_coords_free(&uv);
    }

    /* Flip detection with deliberately mirrored UVs. */
    {
        geomlib_uv_coords uv = geomlib_uv_coords_create(disk.vertex_count, disk.triangle_count);
        for (size_t i = 0; i < uv.vertex_count; ++i) {
            uv.uvs[i].x = (float)i / uv.vertex_count;
            uv.uvs[i].y = (i & 1) ? 1.0f : 0.0f;
        }
        geomlib_param_stats stats;
        geomlib_compute_param_stats(&disk, &uv, &stats);
        check(stats.has_flips, "flip detection works");
        geomlib_uv_coords_free(&uv);
    }

    geomlib_trimesh_free(&disk);
    if (failures != 0) {
        std::printf("%d parameterization test(s) failed\n", failures);
        return 1;
    }
    std::printf("All parameterization tests passed\n");
    return 0;
}
