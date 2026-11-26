/*
 * geomlib - Alpha wrapping (see geomlib_alpha_wrap.h)
 *
 * Clean-room implementation of the alpha-wrap idea from CGAL: build a dense
 * signed-distance field from the input, offset it outward by `alpha`,
 * optionally smooth the field, and extract the watertight isosurface.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_alpha_wrap.h"
#include "geomlib/geomlib_sdf.h"

#include <math.h>
#include <stdlib.h>

geomlib_alpha_wrap_params geomlib_alpha_wrap_params_default(void)
{
    geomlib_alpha_wrap_params p;
    p.resolution = 48;
    p.alpha = 0.0f;
    p.smoothing_steps = 1;
    return p;
}

int geomlib_alpha_wrap(const geomlib_trimesh *input,
                       const geomlib_alpha_wrap_params *params,
                       geomlib_trimesh *out)
{
    geomlib_sdf sdf;
    geomlib_aabb box;
    geomlib_vec3 size;
    float diag, alpha;
    uint32_t resolution;
    int steps, i;
    int rc;

    if (!input || !out || !input->vertices || !input->indices ||
        input->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    box = geomlib_trimesh_compute_aabb(input);
    size = geomlib_aabb_size(box);
    diag = geomlib_vec3_length(size);
    if (!(diag > 0.0f)) return GEOMLIB_MESH_ERROR_INVALID;

    resolution = params ? params->resolution : 48;
    if (resolution < 16) resolution = 16;
    if (resolution > 160) resolution = 160;
    alpha = params ? params->alpha : 0.0f;
    if (!(alpha > 0.0f)) alpha = 0.01f * diag;
    steps = params ? params->smoothing_steps : 1;
    if (steps < 0) steps = 0;
    if (steps > 8) steps = 8;

    {
        geomlib_sdf_params sp = geomlib_sdf_params_default();
        sp.resolution = resolution;
        sp.voxel_size = diag / resolution;
        /* Open inputs have no well-defined inside/outside; an unsigned field
         * offset outward still yields a watertight enclosing slab. */
        sp.compute_sign = false;
        /* Enough padding keeps the offset slab strictly inside the grid so
         * the extracted surface is closed rather than clipped at the border. */
        sp.padding = 8.0f;
        sp.use_bvh = true;
        sp.algorithm = GEOMLIB_SDF_ALGORITHM_BRUTE_FORCE;
        sdf = geomlib_mesh_to_sdf(input, &sp);
    }
    if (!sdf.values) return GEOMLIB_MESH_ERROR_INVALID;

    if (geomlib_sdf_dilate(&sdf, alpha) != GEOMLIB_MESH_SUCCESS) {
        geomlib_sdf_free(&sdf);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    for (i = 0; i < steps; ++i) {
        if (geomlib_sdf_smooth(&sdf, 1.0f) != GEOMLIB_MESH_SUCCESS) break;
    }
    rc = geomlib_sdf_extract_surface(&sdf, out);
    geomlib_sdf_free(&sdf);
    return rc;
}
