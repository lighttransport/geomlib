/*
 * SDF ray casting: rays from outside and inside hit the sphere surface at
 * the expected distance with outward normals; a miss returns hit=0.
 */
#include <geomlib/geomlib_sdf_raycast.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh sphere = {};
    int rc = geomlib_create_sphere(&sphere, 1.0f, 24, 16);
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere created");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_sdf_params p = geomlib_sdf_params_default();
    p.resolution = 48;
    p.compute_sign = true;
    p.use_bvh = true;
    p.padding = 24.0f;
    geomlib_sdf sdf = geomlib_mesh_to_sdf(&sphere, &p);
    check(sdf.values != nullptr, "SDF built");

    geomlib_sdf_ray_hit hit;
    geomlib_ray ray = geomlib_ray_create(geomlib_vec3_create(1.9f,0,0),
                                         geomlib_vec3_create(-1,0,0));
    rc = geomlib_sdf_raycast(&sdf, ray, 10.0f, 256, &hit);
    std::printf("ray outside rc=%d hit=%d t=%.4f radius=%.4f nx=%.3f steps=%u\n",
                rc, hit.hit, hit.t, geomlib_vec3_length(hit.position),
                hit.normal.x, hit.steps);
    check(rc == GEOMLIB_MESH_SUCCESS && hit.hit, "outside ray hits sphere");
    if (hit.hit) {
        check(std::fabs(hit.t-0.9f) < 0.06f, "hit distance is about 0.9");
        check(std::fabs(geomlib_vec3_length(hit.position)-1.0f) < 0.04f,
              "hit point lies on the sphere");
        check(hit.normal.x > 0.9f, "outward normal points along +x");
    }

    ray = geomlib_ray_create(geomlib_vec3_create(0,0,0),
                             geomlib_vec3_create(1,0,0));
    rc = geomlib_sdf_raycast(&sdf, ray, 10.0f, 256, &hit);
    std::printf("ray inside  rc=%d hit=%d t=%.4f\n", rc, hit.hit, hit.t);
    check(rc == GEOMLIB_MESH_SUCCESS && hit.hit && hit.t > 0.9f && hit.t < 1.1f,
          "inside ray exits at radius 1");

    ray = geomlib_ray_create(geomlib_vec3_create(2,2,2),
                             geomlib_vec3_create(1,0,0));
    rc = geomlib_sdf_raycast(&sdf, ray, 10.0f, 256, &hit);
    std::printf("ray miss    rc=%d hit=%d\n", rc, hit.hit);
    check(rc == GEOMLIB_MESH_SUCCESS && !hit.hit, "off-surface ray misses");

    geomlib_sdf_free(&sdf);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d SDF raycast test(s) failed\n", failures);
        return 1;
    }
    std::printf("All SDF raycast tests passed\n");
    return 0;
}
