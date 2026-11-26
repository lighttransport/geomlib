/*
 * RANSAC shape detection tests: two parallel planes in a point cloud.
 */
#include <geomlib/geomlib_shape_detection.h>

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
    for (int i = 0; i < 100; ++i) {
        float x = (float)(i % 10) * 0.1f;
        float y = (float)(i / 10) * 0.1f;
        geomlib_vec3 p1 = {x, y, 0.0f};
        geomlib_vec3 p2 = {x, y, 1.0f};
        geomlib_vec3 n = {0, 0, 1};
        geomlib_pointcloud_add_point(&cloud, p1, &n, 1.0f);
        geomlib_pointcloud_add_point(&cloud, p2, &n, 1.0f);
    }

    size_t count = 0;
    geomlib_detected_shape *shapes = geomlib_detect_shapes_ransac(
        &cloud, GEOMLIB_SHAPE_PLANE, 0.02f, 2, 42, &count);
    std::printf("detected %zu planes\n", count);
    check(shapes != NULL && count == 2, "two planes detected");
    if (shapes && count >= 1) {
        check(shapes[0].inlier_count >= 80, "first plane has most points as inliers");
        check(std::fabs(std::fabs(shapes[0].normal.z) - 1.0f) < 0.05f,
              "detected plane normal is vertical");
    }
    if (shapes && count >= 2) {
        check(shapes[1].inlier_count >= 80, "second plane has most points as inliers");
        check(std::fabs(std::fabs(shapes[1].center.z - shapes[0].center.z) - 1.0f) < 0.05f,
              "plane separation is one");
    }
    geomlib_detected_shapes_free(shapes, count);
    geomlib_pointcloud_free(&cloud);

    if (failures != 0) {
        std::printf("%d shape detection test(s) failed\n", failures);
        return 1;
    }
    std::printf("All shape detection tests passed\n");
    return 0;
}
