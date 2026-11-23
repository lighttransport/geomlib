/* Compare SDF generation modes and zero-level extraction against the dense baseline. */
#include "geomlib/geomlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { repetitions = 3 };

enum sdf_mode {
    MODE_BRUTE_FORCE,
    MODE_FAST_SWEEPING,
    MODE_SPARSE_BAND,
    MODE_FEATURE_EXTRACTION
};

static int compare_time(const void *left, const void *right) {
    double a = *(const double*)left, b = *(const double*)right;
    return (a > b) - (a < b);
}

static const char *mode_name(enum sdf_mode mode) {
    switch (mode) {
        case MODE_BRUTE_FORCE: return "brute_force";
        case MODE_FAST_SWEEPING: return "fast_sweeping";
        case MODE_SPARSE_BAND: return "sparse_band";
        case MODE_FEATURE_EXTRACTION: return "feature_extraction";
    }
    return "unknown";
}

static int benchmark(const char *fixture, const geomlib_trimesh *source, uint32_t resolution,
                     enum sdf_mode mode) {
    double generation[repetitions], extraction[repetitions];
    size_t voxels = 0, stored = 0, bytes = 0, faces = 0;
    for (int i = 0; i < repetitions; i++) {
        geomlib_sdf_params params = geomlib_sdf_params_default();
        params.resolution = resolution;
        params.voxel_size = 0.0f;
        if (mode == MODE_FAST_SWEEPING || mode == MODE_FEATURE_EXTRACTION)
            params.algorithm = GEOMLIB_SDF_ALGORITHM_FAST_SWEEPING;
        if (mode == MODE_SPARSE_BAND) {
            params.narrow_band = true;
            params.band_width = 3.0f;
            params.sparse = true;
        }

        clock_t start = clock();
        geomlib_sdf sdf = geomlib_mesh_to_sdf(source, &params);
        clock_t generated = clock();
        if (!sdf.values && !sdf.band) return 1;

        geomlib_trimesh mesh = {0};
        int status;
        if (mode == MODE_FEATURE_EXTRACTION)
            status = geomlib_sdf_extract_surface_features(&sdf, source, 45.0f, true, &mesh);
        else
            status = geomlib_sdf_extract_surface(&sdf, &mesh);
        clock_t extracted = clock();

        voxels = sdf.total_voxels;
        stored = geomlib_sdf_stored_voxels(&sdf);
        geomlib_sdf_stats stats = geomlib_compute_sdf_stats(&sdf);
        bytes = sdf.is_sparse ? stats.memory_bytes : sdf.total_voxels * sizeof(float);
        faces = mesh.triangle_count;
        geomlib_trimesh_free(&mesh);
        geomlib_sdf_free(&sdf);
        if (status != GEOMLIB_MESH_SUCCESS || !faces || start == (clock_t)-1 ||
            generated == (clock_t)-1 || extracted == (clock_t)-1) return 1;

        generation[i] = 1000.0 * (generated - start) / CLOCKS_PER_SEC;
        extraction[i] = 1000.0 * (extracted - generated) / CLOCKS_PER_SEC;
    }
    qsort(generation, repetitions, sizeof(double), compare_time);
    qsort(extraction, repetitions, sizeof(double), compare_time);
    printf("%s,%s,%u,%zu,%zu,%zu,%.3f,%.3f,%zu\n", fixture, mode_name(mode),
           (unsigned)resolution, voxels, stored, bytes, generation[repetitions / 2],
           extraction[repetitions / 2], faces);
    return 0;
}

int main(void) {
    geomlib_trimesh sphere = {0};
    geomlib_trimesh box = {0};
    if (geomlib_create_sphere(&sphere, 1.0f, 24, 12) != GEOMLIB_MESH_SUCCESS) return 1;
    if (geomlib_create_box(&box, 1.0f, 1.0f, 1.0f) != GEOMLIB_MESH_SUCCESS) return 1;
    printf("Sphere fixture: %zu input triangles\n", sphere.triangle_count);
    printf("Box fixture: %zu input triangles\n", box.triangle_count);
    puts("fixture,mode,resolution,voxels,stored_voxels,grid_bytes,generation_ms,extraction_ms,triangles");

    const uint32_t resolutions[] = {24, 48};
    int failed = 0;
    for (size_t i = 0; i < sizeof(resolutions) / sizeof(resolutions[0]); i++) {
        uint32_t resolution = resolutions[i];
        failed |= benchmark("sphere", &sphere, resolution, MODE_BRUTE_FORCE);
        failed |= benchmark("sphere", &sphere, resolution, MODE_FAST_SWEEPING);
        failed |= benchmark("sphere", &sphere, resolution, MODE_SPARSE_BAND);
        failed |= benchmark("box", &box, resolution, MODE_BRUTE_FORCE);
        failed |= benchmark("box", &box, resolution, MODE_FAST_SWEEPING);
        failed |= benchmark("box", &box, resolution, MODE_FEATURE_EXTRACTION);
    }

    geomlib_trimesh_free(&sphere);
    geomlib_trimesh_free(&box);
    return failed;
}
