/*
 * Mean-value cage deformation tests: affine precision, weight normalization
 * and exact reproduction of deformed cage vertices.
 */
#include <geomlib/geomlib_meancoords.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    std::vector<geomlib_vec2> cage = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1)
    };
    std::vector<geomlib_vec2> q = {
        geomlib_vec2_create(0.2f,0.3f), geomlib_vec2_create(-0.5f,0.1f),
        cage[0], cage[2]
    };
    std::vector<geomlib_vec2> out(q.size());

    int rc = geomlib_deform_mean_value2d(cage.data(), cage.data(), cage.size(),
                                         q.data(), q.size(), out.data());
    std::printf("mvc identity rc=%d first=(%.4f,%.4f)\n", rc, out[0].x, out[0].y);
    check(rc == GEOMLIB_MESH_SUCCESS, "identity deformation succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        double err=0.0;
        for (size_t i=0;i<q.size();++i)
            err += std::hypot((double)out[i].x-q[i].x,(double)out[i].y-q[i].y);
        check(err < 1e-3, "identity cage reproduces query points");
    }

    std::vector<geomlib_vec2> moved(cage.size());
    for (size_t i=0;i<cage.size();++i)
        moved[i]=geomlib_vec2_create(cage[i].x+2.0f,cage[i].y+1.0f);
    rc = geomlib_deform_mean_value2d(cage.data(), moved.data(), cage.size(),
                                     q.data(), q.size(), out.data());
    std::printf("mvc translation first=(%.4f,%.4f)\n", out[0].x, out[0].y);
    check(rc == GEOMLIB_MESH_SUCCESS &&
          std::hypot((double)out[0].x-(q[0].x+2.0),
                     (double)out[0].y-(q[0].y+1.0)) < 1e-3,
          "translation is reproduced exactly");
    check(std::hypot((double)out[2].x-moved[0].x,
                     (double)out[2].y-moved[0].y) < 1e-3,
          "cage vertex maps exactly");

    for (size_t i=0;i<cage.size();++i)
        moved[i]=geomlib_vec2_create(cage[i].x*2.0f,cage[i].y*2.0f);
    rc = geomlib_deform_mean_value2d(cage.data(), moved.data(), cage.size(),
                                     q.data(), q.size(), out.data());
    std::printf("mvc scale first=(%.4f,%.4f)\n", out[0].x, out[0].y);
    check(rc == GEOMLIB_MESH_SUCCESS &&
          std::hypot((double)out[0].x-2.0*q[0].x,
                     (double)out[0].y-2.0*q[0].y) < 1e-3,
          "uniform scale is reproduced exactly");

    double w[4];
    rc = geomlib_mean_value_weights2d(cage.data(), cage.size(), q[0], w);
    double sum=0.0;
    bool positive=true;
    for (int i=0;i<4;++i) { sum+=w[i]; if (w[i] < 0.0) positive=false; }
    std::printf("mvc weights sum=%.6f min=%.4f\n", sum, w[3]);
    check(rc == GEOMLIB_MESH_SUCCESS && std::fabs(sum-1.0) < 1e-6 && positive,
          "weights sum to one and are positive");

    if (failures != 0) {
        std::printf("%d mean value test(s) failed\n", failures);
        return 1;
    }
    std::printf("All mean value tests passed\n");
    return 0;
}
