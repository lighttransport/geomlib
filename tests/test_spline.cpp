/*
 * Catmull-Rom resampling tests: interpolation, collinearity, closed-loop
 * closure and finite output for both parameterizations.
 */
#include <geomlib/geomlib_spline.h>

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
    std::vector<geomlib_vec3> line = {
        geomlib_vec3_create(0,0,0), geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(2,0,0), geomlib_vec3_create(3,0,0)
    };
    size_t count=geomlib_catmull_output_count(line.size(),0,8);
    std::vector<geomlib_vec3> out(count);
    int rc=geomlib_resample_catmull_rom(line.data(),line.size(),0,8,
                                        GEOMLIB_CATMULL_UNIFORM,out.data());
    std::printf("spline open rc=%d count=%zu\n", rc, count);
    check(rc==GEOMLIB_MESH_SUCCESS && count==25, "open output count is 25");
    if (rc==GEOMLIB_MESH_SUCCESS) {
        float maxoff=0.0f;
        for (auto &p:out) maxoff=std::fmax(maxoff,std::fabs(p.y)+std::fabs(p.z));
        check(maxoff<1e-4f,"collinear controls produce a straight spline");
        check(std::fabs(out[0].x)<1e-5f && std::fabs(out.back().x-3.0f)<1e-5f,
              "endpoints are interpolated");
        check(std::fabs(out[8].x-1.0f)<1e-4f,"control point is interpolated");
    }

    std::vector<geomlib_vec3> square = {
        geomlib_vec3_create(1,0,0), geomlib_vec3_create(0,1,0),
        geomlib_vec3_create(-1,0,0), geomlib_vec3_create(0,-1,0)
    };
    count=geomlib_catmull_output_count(square.size(),1,6);
    out.assign(count,geomlib_vec3_create(0,0,0));
    rc=geomlib_resample_catmull_rom(square.data(),square.size(),1,6,
                                    GEOMLIB_CATMULL_CENTRIPETAL,out.data());
    std::printf("spline closed rc=%d count=%zu\n", rc, count);
    check(rc==GEOMLIB_MESH_SUCCESS && count==25,"closed output count is 25");
    if (rc==GEOMLIB_MESH_SUCCESS) {
        bool finite=true;
        for (auto &p:out)
            if (!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)) finite=false;
        check(finite,"all samples finite");
        check(std::fabs(out.front().x-out.back().x)<1e-5f &&
              std::fabs(out.front().y-out.back().y)<1e-5f,"closed loop closes");
        float maxr=0.0f,minr=1e9f;
        for (size_t i=0;i+1<count;++i) {
            float r=geomlib_vec3_length(out[i]);
            maxr=std::fmax(maxr,r); minr=std::fmin(minr,r);
        }
        std::printf("spline radii [%.3f, %.3f]\n",minr,maxr);
        check(maxr<1.3f && minr>0.7f,"samples stay near the circle");
    }

    std::vector<geomlib_vec3> uneven = {
        geomlib_vec3_create(0,0,0), geomlib_vec3_create(0.2f,0,0),
        geomlib_vec3_create(3,0,0), geomlib_vec3_create(3.2f,1,0)
    };
    std::vector<geomlib_vec3> a(geomlib_catmull_output_count(uneven.size(),0,10));
    std::vector<geomlib_vec3> b(a.size());
    geomlib_resample_catmull_rom(uneven.data(),uneven.size(),0,10,
                                 GEOMLIB_CATMULL_UNIFORM,a.data());
    geomlib_resample_catmull_rom(uneven.data(),uneven.size(),0,10,
                                 GEOMLIB_CATMULL_CENTRIPETAL,b.data());
    double diff=0.0;
    for (size_t i=0;i<a.size();++i)
        diff+=geomlib_vec3_length(geomlib_vec3_sub(a[i],b[i]));
    std::printf("spline uniform vs centripetal diff=%.4f\n",diff);
    check(diff>0.01,"parameterizations differ for uneven spacing");

    if (failures!=0) {
        std::printf("%d spline test(s) failed\n",failures);
        return 1;
    }
    std::printf("All spline tests passed\n");
    return 0;
}
