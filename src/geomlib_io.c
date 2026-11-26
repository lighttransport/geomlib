/*
 * geomlib - Mesh file I/O (see geomlib_io.h)
 *
 * Clean-room readers/writers for OBJ, OFF, ASCII PLY and STL. C11 only;
 * no FILE*, no locale-dependent number parsing beyond strtod, wasm32 safe.
 */

#include "geomlib/geomlib_io.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Small dynamic string buffer used by the writers                           */
/* ------------------------------------------------------------------------- */

/* Little-endian 32-bit word from a little-endian byte stream. */
static uint32_t ntohl32(uint32_t v)
{
    unsigned char b[4];
    memcpy(b, &v, 4);
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static uint32_t htonl32(uint32_t v)
{
    return ntohl32(v);
}

typedef struct {
    char *data;
    size_t len;
    size_t cap;
} strbuf;

static int sb_reserve(strbuf *b, size_t extra)
{
    size_t need = b->len + extra + 1;
    char *p;
    if (need <= b->cap) return GEOMLIB_IO_SUCCESS;
    if (b->cap == 0) b->cap = 4096;
    while (b->cap < need) b->cap *= 2;
    p = (char *)realloc(b->data, b->cap);
    if (!p) return GEOMLIB_IO_ERROR_MEMORY;
    b->data = p;
    return GEOMLIB_IO_SUCCESS;
}

static int sb_append(strbuf *b, const char *s)
{
    size_t n = strlen(s);
    if (sb_reserve(b, n) != GEOMLIB_IO_SUCCESS) return GEOMLIB_IO_ERROR_MEMORY;
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
    return GEOMLIB_IO_SUCCESS;
}

static int sb_printf(strbuf *b, const char *fmt, ...)
{
    char tmp[256];
    va_list ap;
    int n;
    va_start(ap, fmt);
    n = vsnprintf(tmp, sizeof(tmp), fmt, ap);
    va_end(ap);
    if (n < 0) return GEOMLIB_IO_ERROR_INVALID;
    if ((size_t)n >= sizeof(tmp)) return GEOMLIB_IO_ERROR_INVALID;
    return sb_append(b, tmp);
}

/* ------------------------------------------------------------------------- */
/* Format detection                                                          */
/* ------------------------------------------------------------------------- */

const char *geomlib_mesh_detect_format(const void *data, size_t size)
{
    const unsigned char *p = (const unsigned char *)data;
    if (!p || size < 4) return NULL;
    if (p[0] == 'p' && p[1] == 'l' && p[2] == 'y' && (p[3] == '\n' || p[3] == '\r')) {
        /* Binary little-endian PLY declares its format in the header. */
        const char *s = (const char *)p;
        const char *e = (const char *)p + size;
        const char *fmt = strstr(s, "format binary_little_endian");
        if (fmt && fmt < e) return "plyb";
        return "ply";
    }
    if (p[0] == 'O' && p[1] == 'F' && p[2] == 'F') return "off";
    if (size >= 84) {
        /* Binary STL: a uint32 triangle count at byte 80 and exactly
         * 50 bytes per triangle follow, independent of the header text. */
        unsigned char hdr[4];
        uint32_t count;
        memcpy(hdr, p + 80, 4);
        count = (uint32_t)hdr[0] | ((uint32_t)hdr[1] << 8) |
                ((uint32_t)hdr[2] << 16) | ((uint32_t)hdr[3] << 24);
        if (count < 1000000 && size == 84 + 50 * (size_t)count) return "stl";
    }
    if (size >= 5 && memcmp(p, "solid", 5) == 0) {
        return "stl"; /* ASCII STL */
    }
    if (size >= 5 && memcmp(p, "solid", 5) == 0) return "stl";
    if (p[0] == 'v' || p[0] == 'f' || p[0] == '#') return "obj";
    return NULL;
}

/* ------------------------------------------------------------------------- */
/* Text helpers                                                              */
/* ------------------------------------------------------------------------- */

typedef struct {
    const char *cur;
    const char *end;
} linescan;

static int next_line(linescan *ls, const char **begin, const char **end)
{
    const char *s = ls->cur;
    const char *e;
    if (s >= ls->end) return 0;
    e = s;
    while (e < ls->end && *e != '\n') ++e;
    *begin = s;
    *end = e;
    ls->cur = (e < ls->end) ? e + 1 : e;
    return 1;
}

static size_t split_tokens(const char *begin, const char *end,
                           const char *tokens[], size_t max_tokens)
{
    size_t n = 0;
    const char *p = begin;
    while (p < end && n < max_tokens) {
        while (p < end && isspace((unsigned char)*p)) ++p;
        if (p >= end) break;
        tokens[n++] = p;
        while (p < end && !isspace((unsigned char)*p)) ++p;
    }
    return n;
}

/* Tokens point into the line buffer and are not NUL-terminated. */
static int token_is(const char *tok, const char *line_end, const char *word)
{
    const char *p = tok;
    size_t i = 0;
    while (p < line_end && !isspace((unsigned char)*p)) { ++p; ++i; }
    return i == strlen(word) && strncmp(tok, word, i) == 0;
}

static double token_double(const char *tok, const char *end)
{
    char tmp[64];
    size_t n = 0;
    while (tok < end && n + 1 < sizeof(tmp) && !isspace((unsigned char)*tok)) {
        tmp[n++] = *tok++;
    }
    tmp[n] = '\0';
    return strtod(tmp, NULL);
}

static long token_index(const char *tok, const char *line_end, long vertex_count)
{
    /* OBJ indices may be "i", "i/vt/vn", "i//vn" and may be negative. */
    char tmp[64];
    size_t n = 0;
    long value;
    while (tok < line_end && n + 1 < sizeof(tmp) &&
           *tok != '/' && !isspace((unsigned char)*tok)) {
        tmp[n++] = *tok++;
    }
    tmp[n] = '\0';
    if (n == 0) return 0;
    value = strtol(tmp, NULL, 10);
    if (value < 0) value = vertex_count + value + 1;
    return value;
}

typedef struct {
    geomlib_vec3 *v;
    size_t count, cap;
} vertex_list;

static int vl_push(vertex_list *l, geomlib_vec3 p)
{
    if (l->count == l->cap) {
        size_t nc = l->cap ? l->cap * 2 : 256;
        geomlib_vec3 *nv = (geomlib_vec3 *)realloc(l->v, nc * sizeof(geomlib_vec3));
        if (!nv) return GEOMLIB_IO_ERROR_MEMORY;
        l->v = nv;
        l->cap = nc;
    }
    l->v[l->count++] = p;
    return GEOMLIB_IO_SUCCESS;
}

static int append_triangle(uint32_t **idx, size_t *icount, size_t *icap,
                           uint32_t a, uint32_t b, uint32_t c)
{
    if (*icount + 3 > *icap) {
        size_t nc = *icap ? *icap * 2 : 256;
        uint32_t *ni = (uint32_t *)realloc(*idx, nc * sizeof(uint32_t));
        if (!ni) return GEOMLIB_IO_ERROR_MEMORY;
        *idx = ni;
        *icap = nc;
    }
    (*idx)[(*icount)++] = a;
    (*idx)[(*icount)++] = b;
    (*idx)[(*icount)++] = c;
    return GEOMLIB_IO_SUCCESS;
}

static int build_from_lists(geomlib_trimesh *mesh, const vertex_list *vl,
                            const uint32_t *idx, size_t icount)
{
    size_t i;
    int rc = geomlib_trimesh_init(mesh, vl->count, icount / 3);
    if (rc != GEOMLIB_MESH_SUCCESS) return GEOMLIB_IO_ERROR_INVALID;
    for (i = 0; i < vl->count; ++i) {
        if (geomlib_trimesh_add_vertex(mesh, vl->v[i], NULL, NULL) < 0) {
            geomlib_trimesh_free(mesh);
            return GEOMLIB_IO_ERROR_INVALID;
        }
    }
    for (i = 0; i + 2 < icount; i += 3) {
        if (geomlib_trimesh_add_triangle(mesh, idx[i], idx[i + 1], idx[i + 2]) < 0) {
            geomlib_trimesh_free(mesh);
            return GEOMLIB_IO_ERROR_INVALID;
        }
    }
    return GEOMLIB_IO_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/* Readers                                                                   */
/* ------------------------------------------------------------------------- */

static int load_obj(geomlib_trimesh *mesh, const char *text, size_t size)
{
    linescan ls = { text, text + size };
    const char *b, *e;
    vertex_list vl = {0};
    uint32_t *idx = NULL;
    size_t icount = 0, icap = 0;
    int rc = GEOMLIB_IO_SUCCESS;
    const char *tok[64];

    while (next_line(&ls, &b, &e)) {
        size_t nt = split_tokens(b, e, tok, 64);
        if (nt == 0 || tok[0] >= e) continue;
        if (token_is(tok[0], e, "v")) {
            geomlib_vec3 p;
            if (nt < 4) { rc = GEOMLIB_IO_ERROR_INVALID; break; }
            p.x = (float)token_double(tok[1], e);
            p.y = (float)token_double(tok[2], e);
            p.z = (float)token_double(tok[3], e);
            rc = vl_push(&vl, p);
            if (rc != GEOMLIB_IO_SUCCESS) break;
        } else if (token_is(tok[0], e, "f")) {
            uint32_t facev[64];
            size_t fn = 0, k;
            if (nt < 4) { rc = GEOMLIB_IO_ERROR_INVALID; break; }
            for (k = 1; k < nt; ++k) {
                long vi = token_index(tok[k], e, (long)vl.count);
                if (vi < 1 || (size_t)vi > vl.count) {
                    rc = GEOMLIB_IO_ERROR_INVALID;
                    break;
                }
                facev[fn++] = (uint32_t)(vi - 1); /* OBJ is 1-based */
            }
            if (rc != GEOMLIB_IO_SUCCESS) break;
            for (k = 1; k + 1 < fn; ++k) {
                rc = append_triangle(&idx, &icount, &icap, facev[0], facev[k], facev[k + 1]);
                if (rc != GEOMLIB_IO_SUCCESS) break;
            }
            if (rc != GEOMLIB_IO_SUCCESS) break;
        }
    }

    if (rc == GEOMLIB_IO_SUCCESS) {
        rc = build_from_lists(mesh, &vl, idx, icount);
    }
    free(vl.v);
    free(idx);
    return rc;
}

static int load_off(geomlib_trimesh *mesh, const char *text, size_t size)
{
    linescan ls = { text, text + size };
    const char *b, *e;
    vertex_list vl = {0};
    uint32_t *idx = NULL;
    size_t icount = 0, icap = 0;
    size_t nv = 0, nf = 0, i;
    int rc = GEOMLIB_IO_SUCCESS;
    const char *tok[64];
    int header_seen = 0;

    while (next_line(&ls, &b, &e)) {
        size_t nt = split_tokens(b, e, tok, 64);
        if (nt == 0) continue;
        if (!header_seen) {
            if (tok[0][0] != 'O' || tok[0][1] != 'F' || tok[0][2] != 'F') {
                return GEOMLIB_IO_ERROR_INVALID;
            }
            header_seen = 1;
            continue;
        }
        if (nv == 0 && nf == 0) {
            if (nt < 2) return GEOMLIB_IO_ERROR_INVALID;
            nv = (size_t)strtoul(tok[0], NULL, 10);
            nf = (size_t)strtoul(tok[1], NULL, 10);
            continue;
        }
        if (vl.count < nv) {
            geomlib_vec3 p;
            if (nt < 3) { rc = GEOMLIB_IO_ERROR_INVALID; break; }
            p.x = (float)token_double(tok[0], e);
            p.y = (float)token_double(tok[1], e);
            p.z = (float)token_double(tok[2], e);
            rc = vl_push(&vl, p);
            if (rc != GEOMLIB_IO_SUCCESS) break;
        } else {
            size_t k, count;
            if (nt == 0) { rc = GEOMLIB_IO_ERROR_INVALID; break; }
            count = (size_t)strtoul(tok[0], NULL, 10);
            if (count < 3 || count > 63 || nt < count + 1) {
                rc = GEOMLIB_IO_ERROR_INVALID;
                break;
            }
            {
                uint32_t first = (uint32_t)strtoul(tok[1], NULL, 10);
                for (k = 2; k < count; ++k) {
                    uint32_t a = (uint32_t)strtoul(tok[k], NULL, 10);
                    uint32_t c = (uint32_t)strtoul(tok[k + 1], NULL, 10);
                    rc = append_triangle(&idx, &icount, &icap, first, a, c);
                    if (rc != GEOMLIB_IO_SUCCESS) break;
                }
            }
            if (rc != GEOMLIB_IO_SUCCESS) break;
        }
    }
    (void)i;
    if (rc == GEOMLIB_IO_SUCCESS) rc = build_from_lists(mesh, &vl, idx, icount);
    free(vl.v);
    free(idx);
    return rc;
}

static int load_ply(geomlib_trimesh *mesh, const char *text, size_t size)
{
    linescan ls = { text, text + size };
    const char *b, *e;
    const char *tok[64];
    size_t nv = 0, nf = 0, i;
    vertex_list vl = {0};
    uint32_t *idx = NULL;
    size_t icount = 0, icap = 0;
    int rc = GEOMLIB_IO_SUCCESS;
    int header = 1;

    while (next_line(&ls, &b, &e)) {
        size_t nt = split_tokens(b, e, tok, 64);
        if (nt == 0) continue;
        if (header) {
            if (token_is(tok[0], e, "ply")) continue;
            if (token_is(tok[0], e, "format")) {
                if (nt < 2 || !token_is(tok[1], e, "ascii"))
                    return GEOMLIB_IO_ERROR_UNSUPPORTED;
            } else if (token_is(tok[0], e, "element") && nt >= 3) {
                if (token_is(tok[1], e, "vertex")) nv = (size_t)strtoul(tok[2], NULL, 10);
                else if (token_is(tok[1], e, "face")) nf = (size_t)strtoul(tok[2], NULL, 10);
            } else if (token_is(tok[0], e, "end_header")) {
                header = 0;
            }
            continue;
        }

        if (vl.count < nv) {
            geomlib_vec3 p;
            if (nt < 3) { rc = GEOMLIB_IO_ERROR_INVALID; break; }
            p.x = (float)token_double(tok[0], e);
            p.y = (float)token_double(tok[1], e);
            p.z = (float)token_double(tok[2], e);
            rc = vl_push(&vl, p);
            if (rc != GEOMLIB_IO_SUCCESS) break;
        } else {
            size_t k, count;
            if (nt < 4) { rc = GEOMLIB_IO_ERROR_INVALID; break; }
            count = (size_t)strtoul(tok[0], NULL, 10);
            if (count < 3 || count > 63 || nt < count + 1) {
                rc = GEOMLIB_IO_ERROR_INVALID;
                break;
            }
            {
                uint32_t first = (uint32_t)strtoul(tok[1], NULL, 10);
                for (k = 2; k < count; ++k) {
                    uint32_t a = (uint32_t)strtoul(tok[k], NULL, 10);
                    uint32_t c = (uint32_t)strtoul(tok[k + 1], NULL, 10);
                    rc = append_triangle(&idx, &icount, &icap, first, a, c);
                    if (rc != GEOMLIB_IO_SUCCESS) break;
                }
            }
            if (rc != GEOMLIB_IO_SUCCESS) break;
        }
    }
    (void)nf; (void)i;
    if (rc == GEOMLIB_IO_SUCCESS) rc = build_from_lists(mesh, &vl, idx, icount);
    free(vl.v);
    free(idx);
    return rc;
}

static int load_stl(geomlib_trimesh *mesh, const void *data, size_t size)
{
    const unsigned char *p = (const unsigned char *)data;
    vertex_list vl = {0};
    uint32_t *idx = NULL;
    size_t icount = 0, icap = 0, tri, i;
    int rc = GEOMLIB_IO_SUCCESS;

    /* Binary STL: 80-byte header, uint32 triangle count, 50 bytes/triangle. */
    if (size >= 84) {
        unsigned char hdr[4];
        uint32_t count;
        memcpy(hdr, p + 80, 4);
        count = (uint32_t)hdr[0] | ((uint32_t)hdr[1] << 8) |
                ((uint32_t)hdr[2] << 16) | ((uint32_t)hdr[3] << 24);
        if (size == 84 + 50 * (size_t)count) {
            rc = geomlib_trimesh_init(mesh, (size_t)count * 3, count);
            if (rc != GEOMLIB_MESH_SUCCESS) return GEOMLIB_IO_ERROR_INVALID;
            for (tri = 0; tri < count; ++tri) {
                const unsigned char *rec = p + 84 + 50 * tri;
                float verts[9];
                uint32_t ids[3];
                for (i = 0; i < 9; ++i) {
                    unsigned char raw[4];
                    memcpy(raw, rec + 12 + i * 4, 4);
                    memcpy(&verts[i], raw, 4);
                }
                for (i = 0; i < 3; ++i) {
                    ids[i] = (uint32_t)geomlib_trimesh_add_vertex(mesh,
                        geomlib_vec3_create(verts[i*3], verts[i*3+1], verts[i*3+2]),
                        NULL, NULL);
                }
                geomlib_trimesh_add_triangle(mesh, ids[0], ids[1], ids[2]);
            }
            return GEOMLIB_IO_SUCCESS;
        }
    }

    /* ASCII STL: collect every "vertex x y z". */
    {
        linescan ls = { (const char *)data, (const char *)data + size };
        const char *b, *e;
        const char *tok[64];
        geomlib_vec3 v[3];
        int nv = 0;
        while (next_line(&ls, &b, &e)) {
            size_t nt = split_tokens(b, e, tok, 64);
            if (nt >= 4 && token_is(tok[0], e, "vertex")) {
                v[nv].x = (float)token_double(tok[1], e);
                v[nv].y = (float)token_double(tok[2], e);
                v[nv].z = (float)token_double(tok[3], e);
                if (++nv == 3) {
                    uint32_t ids[3];
                    for (i = 0; i < 3; ++i) {
                        if ((rc = vl_push(&vl, v[i])) != GEOMLIB_IO_SUCCESS) goto done;
                    }
                    ids[0] = (uint32_t)(vl.count - 3);
                    ids[1] = (uint32_t)(vl.count - 2);
                    ids[2] = (uint32_t)(vl.count - 1);
                    if (icount + 3 > icap) {
                        size_t nc = icap ? icap * 2 : 256;
                        uint32_t *ni = (uint32_t *)realloc(idx, nc * sizeof(uint32_t));
                        if (!ni) { rc = GEOMLIB_IO_ERROR_MEMORY; goto done; }
                        idx = ni; icap = nc;
                    }
                    idx[icount++] = ids[0];
                    idx[icount++] = ids[1];
                    idx[icount++] = ids[2];
                    nv = 0;
                }
            }
        }
    }

done:
    if (rc == GEOMLIB_IO_SUCCESS) rc = build_from_lists(mesh, &vl, idx, icount);
    free(vl.v);
    free(idx);
    return rc;
}

/* Binary little-endian PLY with float32 x/y/z and list-of-uchar face indices. */
static int load_ply_binary(geomlib_trimesh *mesh, const char *text, size_t size)
{
    linescan ls = { text, text + size };
    const char *b, *e, *data_start = NULL;
    size_t nv = 0, nf = 0, i, off;
    vertex_list vl = {0};
    uint32_t *idx = NULL;
    size_t icount = 0, icap = 0;
    const unsigned char *cur;
    const char *tok[64];
    int rc = GEOMLIB_IO_SUCCESS;

    while (next_line(&ls, &b, &e)) {
        size_t nt = split_tokens(b, e, tok, 64);
        if (nt == 0) continue;
        if (token_is(tok[0], e, "ply")) continue;
        if (token_is(tok[0], e, "format")) {
            if (nt < 2 || !token_is(tok[1], e, "binary_little_endian"))
                return GEOMLIB_IO_ERROR_UNSUPPORTED;
        } else if (token_is(tok[0], e, "element") && nt >= 3) {
            if (token_is(tok[1], e, "vertex")) nv = (size_t)strtoul(tok[2], NULL, 10);
            else if (token_is(tok[1], e, "face")) nf = (size_t)strtoul(tok[2], NULL, 10);
        } else if (token_is(tok[0], e, "end_header")) {
            data_start = e + 1;
            break;
        }
    }
    if (!data_start || !nv) return GEOMLIB_IO_ERROR_INVALID;
    off = (size_t)(data_start - text);
    if (off + nv*12 + nf*4 > size) return GEOMLIB_IO_ERROR_INVALID;
    cur = (const unsigned char *)text + off;
    for (i = 0; i < nv; ++i) {
        float v[3];
        memcpy(v, cur + i*12, 12);
        if ((rc = vl_push(&vl, geomlib_vec3_create(v[0],v[1],v[2]))) != GEOMLIB_IO_SUCCESS)
            goto done;
    }
    cur += nv*12;
    for (i = 0; i < nf; ++i) {
        uint32_t n = cur[0];
        const uint32_t *fi;
        size_t k;
        cur++;
        if (n < 3 || n > 63) { rc = GEOMLIB_IO_ERROR_INVALID; goto done; }
        fi = (const uint32_t *)cur;
        for (k = 1; k + 1 < n; ++k) {
            uint32_t a = ntohl32(fi[0]);
            uint32_t c = ntohl32(fi[k]);
            if ((rc = append_triangle(&idx, &icount, &icap, a, c, ntohl32(fi[k+1]))) != GEOMLIB_IO_SUCCESS)
                goto done;
        }
        cur += n*4;
    }

done:
    if (rc == GEOMLIB_IO_SUCCESS) rc = build_from_lists(mesh, &vl, idx, icount);
    free(vl.v);
    free(idx);
    return rc;
}

int geomlib_mesh_load_mem(geomlib_trimesh *mesh, const char *format,
                          const void *data, size_t size)
{
    char *text;
    int rc;
    if (!mesh || !data || size == 0) return GEOMLIB_IO_ERROR_INVALID;
    if (!format) format = geomlib_mesh_detect_format(data, size);
    if (!format) return GEOMLIB_IO_ERROR_INVALID;

    if (strcmp(format, "stl") == 0 && size >= 84) {
        /* Binary STL is parsed directly; ASCII STL uses a text copy. */
        unsigned char hdr[4];
        uint32_t count;
        memcpy(hdr, (const unsigned char *)data + 80, 4);
        count = (uint32_t)hdr[0] | ((uint32_t)hdr[1] << 8) |
                ((uint32_t)hdr[2] << 16) | ((uint32_t)hdr[3] << 24);
        if (size == 84 + 50 * (size_t)count) {
            return load_stl(mesh, data, size);
        }
    }

    text = (char *)malloc(size + 1);
    if (!text) return GEOMLIB_IO_ERROR_MEMORY;
    memcpy(text, data, size);
    text[size] = '\0';
    if (strcmp(format, "obj") == 0) rc = load_obj(mesh, text, size);
    else if (strcmp(format, "off") == 0) rc = load_off(mesh, text, size);
    else if (strcmp(format, "ply") == 0) rc = load_ply(mesh, text, size);
    else if (strcmp(format, "plyb") == 0) rc = load_ply_binary(mesh, text, size);
    else if (strcmp(format, "stl") == 0) rc = load_stl(mesh, text, size);
    else rc = GEOMLIB_IO_ERROR_UNSUPPORTED;
    free(text);
    return rc;
}

/* ------------------------------------------------------------------------- */
/* Writers                                                                   */
/* ------------------------------------------------------------------------- */

int geomlib_mesh_save_mem(const geomlib_trimesh *mesh, const char *format,
                          void **out_data, size_t *out_size)
{
    strbuf b = {0};
    size_t i;
    int rc = GEOMLIB_IO_SUCCESS;
    if (!mesh || !format || !out_data || !out_size) return GEOMLIB_IO_ERROR_INVALID;
    if (!mesh->vertices || !mesh->indices) return GEOMLIB_IO_ERROR_INVALID;

    if (strcmp(format, "obj") == 0) {
        for (i = 0; i < mesh->vertex_count; ++i) {
            geomlib_vec3 v = mesh->vertices[i];
            rc = sb_printf(&b, "v %.9g %.9g %.9g\n", v.x, v.y, v.z);
            if (rc != GEOMLIB_IO_SUCCESS) goto done;
        }
        for (i = 0; i < mesh->triangle_count; ++i) {
            rc = sb_printf(&b, "f %u %u %u\n",
                           mesh->indices[i*3] + 1,
                           mesh->indices[i*3+1] + 1,
                           mesh->indices[i*3+2] + 1);
            if (rc != GEOMLIB_IO_SUCCESS) goto done;
        }
    } else if (strcmp(format, "off") == 0) {
        rc = sb_printf(&b, "OFF\n%zu %zu 0\n", mesh->vertex_count, mesh->triangle_count);
        if (rc != GEOMLIB_IO_SUCCESS) goto done;
        for (i = 0; i < mesh->vertex_count; ++i) {
            geomlib_vec3 v = mesh->vertices[i];
            rc = sb_printf(&b, "%.9g %.9g %.9g\n", v.x, v.y, v.z);
            if (rc != GEOMLIB_IO_SUCCESS) goto done;
        }
        for (i = 0; i < mesh->triangle_count; ++i) {
            rc = sb_printf(&b, "3 %u %u %u\n",
                           mesh->indices[i*3], mesh->indices[i*3+1], mesh->indices[i*3+2]);
            if (rc != GEOMLIB_IO_SUCCESS) goto done;
        }
    } else if (strcmp(format, "ply") == 0) {
        rc = sb_append(&b, "ply\nformat ascii 1.0\n");
        if (rc != GEOMLIB_IO_SUCCESS) goto done;
        rc = sb_printf(&b, "element vertex %zu\n", mesh->vertex_count);
        if (rc != GEOMLIB_IO_SUCCESS) goto done;
        rc = sb_append(&b, "property float x\nproperty float y\nproperty float z\n");
        if (rc != GEOMLIB_IO_SUCCESS) goto done;
        rc = sb_printf(&b, "element face %zu\n", mesh->triangle_count);
        if (rc != GEOMLIB_IO_SUCCESS) goto done;
        rc = sb_append(&b, "property list uchar int vertex_indices\nend_header\n");
        if (rc != GEOMLIB_IO_SUCCESS) goto done;
        for (i = 0; i < mesh->vertex_count; ++i) {
            geomlib_vec3 v = mesh->vertices[i];
            rc = sb_printf(&b, "%.9g %.9g %.9g\n", v.x, v.y, v.z);
            if (rc != GEOMLIB_IO_SUCCESS) goto done;
        }
        for (i = 0; i < mesh->triangle_count; ++i) {
            rc = sb_printf(&b, "3 %u %u %u\n",
                           mesh->indices[i*3], mesh->indices[i*3+1], mesh->indices[i*3+2]);
            if (rc != GEOMLIB_IO_SUCCESS) goto done;
        }
    } else if (strcmp(format, "plyb") == 0) {
        const char header[] =
            "ply\nformat binary_little_endian 1.0\n"
            "element vertex %zu\n"
            "property float x\nproperty float y\nproperty float z\n"
            "element face %zu\n"
            "property list uchar int vertex_indices\nend_header\n";
        size_t hlen, total, p, k;
        unsigned char *bin;
        char htext[512];
        hlen = (size_t)snprintf(htext, sizeof(htext), header,
                                mesh->vertex_count, mesh->triangle_count);
        total = hlen + mesh->vertex_count*12 + mesh->triangle_count*16;
        bin = (unsigned char *)malloc(total);
        if (!bin) { rc = GEOMLIB_IO_ERROR_MEMORY; goto done; }
        memcpy(bin, htext, hlen);
        p = hlen;
        for (i = 0; i < mesh->vertex_count; ++i) {
            memcpy(bin+p, &mesh->vertices[i].x, 4);
            memcpy(bin+p+4, &mesh->vertices[i].y, 4);
            memcpy(bin+p+8, &mesh->vertices[i].z, 4);
            p += 12;
        }
        for (i = 0; i < mesh->triangle_count; ++i) {
            uint32_t n = 3;
            bin[p++] = (unsigned char)n;
            for (k = 0; k < 3; ++k) {
                uint32_t v = htonl32(mesh->indices[i*3+k]);
                memcpy(bin+p, &v, 4);
                p += 4;
            }
        }
        *out_data = bin;
        *out_size = total;
        return GEOMLIB_IO_SUCCESS;
    } else if (strcmp(format, "stl") == 0) {
        size_t total = 84 + 50 * mesh->triangle_count;
        unsigned char *bin = (unsigned char *)calloc(total, 1);
        uint32_t count = (uint32_t)mesh->triangle_count;
        if (!bin) { rc = GEOMLIB_IO_ERROR_MEMORY; goto done; }
        memcpy(bin, "geomlib binary STL", 18);
        bin[80] = (unsigned char)(count & 0xff);
        bin[81] = (unsigned char)((count >> 8) & 0xff);
        bin[82] = (unsigned char)((count >> 16) & 0xff);
        bin[83] = (unsigned char)((count >> 24) & 0xff);
        for (i = 0; i < mesh->triangle_count; ++i) {
            unsigned char *rec = bin + 84 + 50 * i;
            size_t k;
            for (k = 0; k < 3; ++k) {
                geomlib_vec3 v = mesh->vertices[mesh->indices[i*3+k]];
                memcpy(rec + 12 + k*12, &v.x, 4);
                memcpy(rec + 16 + k*12, &v.y, 4);
                memcpy(rec + 20 + k*12, &v.z, 4);
            }
        }
        *out_data = bin;
        *out_size = total;
        return GEOMLIB_IO_SUCCESS;
    } else {
        rc = GEOMLIB_IO_ERROR_UNSUPPORTED;
        goto done;
    }

    *out_data = b.data;
    *out_size = b.len;
    return GEOMLIB_IO_SUCCESS;

done:
    free(b.data);
    return rc;
}
