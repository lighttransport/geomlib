/*
 * 2D polygon Boolean tests over overlapping and contained squares.
 */
#include <geomlib/geomlib_polybool2d.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static double total_area(const geomlib_poly_bool_result &r)
{
    double sum = 0.0;
    for (size_t l = 0; l < r.loop_count; ++l) {
        size_t a=r.offsets[l], b=r.offsets[l+1];
        double s = 0.0;
        for (size_t i = a; i < b; ++i) {
            size_t j = (i+1<b) ? i+1 : a;
            s += (double)r.loops[i].x*r.loops[j].y
               - (double)r.loops[j].x*r.loops[i].y;
        }
        sum += std::fabs(0.5*s);
    }
    return sum;
}

static void run(const geomlib_vec2 *a, size_t na, const geomlib_vec2 *b, size_t nb,
                geomlib_poly_bool_op op, const char *name, double expect)
{
    geomlib_poly_bool_result r = geomlib_poly_bool_result_create();
    int rc = geomlib_polygon_boolean(a,na,b,nb,op,&r);
    double area = rc == GEOMLIB_MESH_SUCCESS ? total_area(r) : 0.0;
    std::printf("polybool %s rc=%d loops=%zu area=%.4f expect=%.4f\n",
                name, rc, r.loop_count, area, expect);
    check(rc == GEOMLIB_MESH_SUCCESS, name);
    if (rc == GEOMLIB_MESH_SUCCESS)
        check(std::fabs(area-expect) < 0.05, name);
    geomlib_poly_bool_result_free(&r);
}

int main(void)
{
    geomlib_vec2 a[4] = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2), geomlib_vec2_create(0,2)
    };
    geomlib_vec2 b[4] = {
        geomlib_vec2_create(1,1), geomlib_vec2_create(3,1),
        geomlib_vec2_create(3,3), geomlib_vec2_create(1,3)
    };
    run(a,4,b,4,GEOMLIB_POLY_UNION,"union of overlapping squares",7.0);
    run(a,4,b,4,GEOMLIB_POLY_INTERSECTION,"intersection of overlapping squares",1.0);
    run(a,4,b,4,GEOMLIB_POLY_DIFFERENCE,"difference of overlapping squares",3.0);

    geomlib_vec2 inner[4] = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1)
    };
    geomlib_vec2 outer[4] = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(2,-1),
        geomlib_vec2_create(2,2), geomlib_vec2_create(-1,2)
    };
    run(inner,4,outer,4,GEOMLIB_POLY_INTERSECTION,"contained intersection",1.0);
    run(inner,4,outer,4,GEOMLIB_POLY_UNION,"contained union",9.0);
    run(inner,4,outer,4,GEOMLIB_POLY_DIFFERENCE,"contained difference",0.0);

    if (failures != 0) {
        std::printf("%d polybool test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polybool tests passed\n");
    return 0;
}
