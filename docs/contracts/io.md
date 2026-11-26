# Contract: mesh I/O

**Status:** implemented (Wave 0). **Sources:** published OBJ, OFF, ASCII PLY
and STL format specifications. No library source consulted.

## Specification

- `geomlib_mesh_load_mem` / `geomlib_mesh_save_mem` operate purely on memory.
- OBJ: vertices, faces with `v/vt/vn` indices and negative indices; n-gon faces
  fan-triangulated. OFF: standard header and face records. PLY: ASCII vertex
  and face elements. STL: ASCII and binary (little-endian).
- Auto-detection from buffer contents.

## Derivation

Format grammars from the public specifications. Binary STL triangle count is
read at byte 80 and validated against the buffer length.

## API / constraints

C11, no `FILE*`, no filesystem dependency — this is the key WASM enabler.
WASM: yes.

Binary little-endian PLY is available under the `"plyb"` format name with the
standard float32 x/y/z and uchar list face properties; auto-detection recognizes
the `format binary_little_endian` header.

## Acceptance

`tests/test_io.cpp`: round-trips five formats (including binary PLY), OBJ
polygon/negative-index case, ASCII STL; `web/smoke.mjs` checks OBJ/PLY and
binary PLY round-trips in Node.

## Lean target

Not applicable (format parsing, not a geometric theorem).
