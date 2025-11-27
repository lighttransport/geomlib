/*
 * geomlib - Texture Atlas Generation
 * Automatic UV atlas generation with chart segmentation and packing
 */

#ifndef GEOMLIB_ATLAS_H
#define GEOMLIB_ATLAS_H

#include "geomlib_mesh.h"
#include "geomlib_parameterization.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Chart segmentation methods */
typedef enum geomlib_atlas_segmentation {
    GEOMLIB_ATLAS_SEG_AUTO,           /* Automatic based on curvature/boundaries */
    GEOMLIB_ATLAS_SEG_PLANARITY,      /* Segment by planarity */
    GEOMLIB_ATLAS_SEG_NORMAL_CLUSTER, /* Segment by normal clustering */
    GEOMLIB_ATLAS_SEG_DISTORTION      /* Segment by parameterization distortion */
} geomlib_atlas_segmentation;

/* Chart packing methods */
typedef enum geomlib_atlas_packing {
    GEOMLIB_ATLAS_PACK_SHELF,         /* Shelf packing (bottom-left) */
    GEOMLIB_ATLAS_PACK_MAXRECTS,      /* MaxRects rectangle packing */
    GEOMLIB_ATLAS_PACK_GUILLOTINE     /* Guillotine packing */
} geomlib_atlas_packing;

/* Chart structure */
typedef struct geomlib_atlas_chart {
    uint32_t *triangles;              /* Triangle indices in this chart */
    size_t triangle_count;
    geomlib_vec2 *uvs;                /* UV coordinates for chart vertices */
    size_t vertex_count;
    uint32_t *vertex_map;             /* Maps chart vertices to mesh vertices */

    /* Bounding box in UV space */
    geomlib_vec2 uv_min;
    geomlib_vec2 uv_max;

    /* Placement in atlas */
    geomlib_vec2 atlas_offset;
    geomlib_vec2 atlas_scale;
} geomlib_atlas_chart;

/* Atlas generation parameters */
typedef struct geomlib_atlas_params {
    geomlib_atlas_segmentation segmentation;
    geomlib_atlas_packing packing;
    geomlib_param_method parameterization;  /* LSCM or Harmonic */

    /* Segmentation parameters */
    float max_chart_area;             /* Maximum relative chart area (0-1) */
    float normal_threshold;           /* Normal angle threshold (degrees) */
    float distortion_threshold;       /* Max stretch ratio for chart */

    /* Packing parameters */
    float gutter_width;               /* Space between charts (in atlas units) */
    int atlas_resolution;             /* Target atlas resolution (power of 2) */
    bool normalize_charts;            /* Normalize charts to similar texel density */

    /* Quality parameters */
    int max_iterations;               /* Max iterations for refinement */
    bool optimize_rotation;           /* Rotate charts for better packing */
} geomlib_atlas_params;

/* Atlas generation result */
typedef struct geomlib_atlas {
    geomlib_atlas_chart *charts;
    size_t chart_count;

    /* Final UVs for entire mesh (stored in mesh->uvs) */
    geomlib_vec2 *uvs;
    size_t uv_count;

    /* Atlas statistics */
    float coverage;                   /* Percentage of atlas used */
    float avg_stretch;                /* Average stretch ratio */
    float max_stretch;                /* Maximum stretch ratio */
} geomlib_atlas;

/* Default atlas parameters */
geomlib_atlas_params geomlib_atlas_params_default(void);

/* Create atlas structure */
geomlib_atlas geomlib_atlas_create(size_t max_charts);

/* Free atlas */
void geomlib_atlas_free(geomlib_atlas *atlas);

/* Generate texture atlas for mesh */
int geomlib_generate_atlas(
    geomlib_trimesh *mesh,
    const geomlib_atlas_params *params,
    geomlib_atlas *atlas
);

/* Chart segmentation - segment mesh into charts */
int geomlib_atlas_segment_mesh(
    const geomlib_trimesh *mesh,
    const geomlib_atlas_params *params,
    geomlib_atlas *atlas
);

/* Parameterize all charts */
int geomlib_atlas_parameterize_charts(
    const geomlib_trimesh *mesh,
    geomlib_atlas *atlas,
    geomlib_param_method method
);

/* Pack charts into atlas space */
int geomlib_atlas_pack_charts(
    geomlib_atlas *atlas,
    const geomlib_atlas_params *params
);

/* Apply atlas UVs to mesh */
int geomlib_atlas_apply_to_mesh(
    geomlib_trimesh *mesh,
    const geomlib_atlas *atlas
);

/* Compute atlas statistics */
void geomlib_atlas_compute_stats(
    const geomlib_trimesh *mesh,
    geomlib_atlas *atlas
);

/* Get chart for triangle */
int geomlib_atlas_get_chart_for_triangle(
    const geomlib_atlas *atlas,
    uint32_t triangle_idx
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ATLAS_H */
