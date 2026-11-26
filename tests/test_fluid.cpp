/*
 * Fluid solver test: pressure projection reduces interior divergence of a
 * rotational velocity field and advection keeps the field finite.
 *
 * The collocated grid cannot cancel divergence exactly in boundary cells,
 * so the property is checked on the interior (as is standard for this
 * discretization).
 */
#include <geomlib/geomlib_fluid.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

/* RMS forward-difference divergence over cells with a full 7-point stencil. */
static float interior_divergence_rms(const geomlib_fluid *f)
{
    double s = 0.0;
    int count = 0;
    float inv = 1.0f/f->dx;
    auto at = [&](const float *a, uint32_t x, uint32_t y, uint32_t z)->size_t {
        return x + f->nx*(y + f->ny*z);
    };
    for (uint32_t z = 1; z+1 < f->nz; ++z)
        for (uint32_t y = 1; y+1 < f->ny; ++y)
            for (uint32_t x = 1; x+1 < f->nx; ++x) {
                float du = f->u[at(f->u,x+1,y,z)] - f->u[at(f->u,x,y,z)];
                float dv = f->v[at(f->v,x,y+1,z)] - f->v[at(f->v,x,y,z)];
                float dw = f->w[at(f->w,x,y,z+1)] - f->w[at(f->w,x,y,z)];
                float d = (du+dv+dw)*inv;
                s += d*d;
                count++;
            }
    return count ? (float)std::sqrt(s/count) : 0.0f;
}

int main(void)
{
    geomlib_fluid f = geomlib_fluid_create(8,8,8,0.1f);
    check(f.u != nullptr, "fluid grid created");
    for (uint32_t z = 0; z < f.nz; ++z)
        for (uint32_t y = 0; y < f.ny; ++y)
            for (uint32_t x = 0; x < f.nx; ++x) {
                size_t id = x + f.nx*(y + f.ny*z);
                f.u[id] = sinf(0.7f*(float)x);
                f.v[id] = -cosf(0.7f*(float)y);
                f.w[id] = 0.0f;
            }
    float before = interior_divergence_rms(&f);
    int rc = geomlib_fluid_project(&f, 1.0f);
    float after = interior_divergence_rms(&f);
    std::printf("fluid project rc=%d interior_div before=%.4f after=%.4f\n",
                rc, before, after);
    check(rc == GEOMLIB_MESH_SUCCESS, "projection succeeds");
    check(after < before*0.1f, "interior divergence reduced by >90%");

    rc = geomlib_fluid_advect(&f, 0.05f);
    bool finite = true;
    for (size_t i = 0; i < (size_t)f.nx*f.ny*f.nz; ++i)
        if (!std::isfinite(f.u[i]) || !std::isfinite(f.v[i]) || !std::isfinite(f.w[i]))
            finite = false;
    std::printf("fluid advect rc=%d finite=%d\n", rc, finite ? 1 : 0);
    check(rc == GEOMLIB_MESH_SUCCESS && finite, "advection keeps velocity finite");
    geomlib_fluid_free(&f);
    if (failures != 0) {
        std::printf("%d fluid test(s) failed\n", failures);
        return 1;
    }
    std::printf("All fluid tests passed\n");
    return 0;
}
