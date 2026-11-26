/*
 * Point-set processing tests: bilateral denoising reduces noise on a plane
 * while outlier removal drops far-away points.
 */
#include <geomlib/geomlib_pointset.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud, 200);
    for (int y = 0; y < 10; ++y)
        for (int x = 0; x < 10; ++x) {
            geomlib_vec3 p = {
                (float)x*0.1f, (float)y*0.1f,
                (((x*7+y*13)%5)-2)*0.01f
            };
            geomlib_vec3 n = {0,0,1};
            geomlib_pointcloud_add_point(&cloud, p, &n, 1.0f);
        }
    double before = 0.0;
    for (size_t i = 0; i < cloud.point_count; ++i)
        before += std::fabs(cloud.points[i].z);
    before /= cloud.point_count;

    geomlib_bilateral_params params = geomlib_bilateral_params_default();
    params.iterations = 2;
    params.neighbors = 12;
    int rc = geomlib_pointcloud_bilateral_smooth(&cloud, &params);
    std::printf("bilateral rc=%d mean|z| %.5f -> %.5f\n",
                rc, before, 0.0);
    check(rc == GEOMLIB_MESH_SUCCESS, "bilateral smoothing succeeds");
    double after = 0.0;
    for (size_t i = 0; i < cloud.point_count; ++i)
        after += std::fabs(cloud.points[i].z);
    after /= cloud.point_count;
    std::printf("bilateral mean|z| %.5f -> %.5f\n", before, after);
    check(after < before, "bilateral smoothing reduces noise");
    check(after < 0.01, "denoised points lie near the plane");

    geomlib_vec3 outlier = {5, 5, 5};
    geomlib_vec3 n = {0,0,1};
    geomlib_pointcloud_add_point(&cloud, outlier, &n, 1.0f);
    size_t removed = geomlib_pointcloud_remove_outliers(&cloud, 0.05f, 8);
    std::printf("outliers removed=%zu\n", removed);
    check(removed >= 1, "outlier removal drops points");
    bool outlier_gone = true;
    for (size_t i = 0; i < cloud.point_count; ++i)
        if (cloud.points[i].x > 2.0f || cloud.points[i].y > 2.0f) outlier_gone = false;
    check(outlier_gone, "the distant outlier was removed");

    geomlib_pointcloud_free(&cloud);
    if (failures != 0) {
        std::printf("%d point-set test(s) failed\n", failures);
        return 1;
    }
    std::printf("All point-set tests passed\n");
    return 0;
}
