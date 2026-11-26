#!/bin/bash
# Check all core-Lean formalization files and verify the generated C table.
set -e
cd "$(dirname "$0")/.."
for f in formal/GeomLibFormal/*.lean; do
    echo "checking $f"
    lean "$f"
done
tmp=$(mktemp)
lean formal/gen_mc_table.lean > "$tmp"
cmp "$tmp" include/geomlib/generated/mc_table.h
rm -f "$tmp"
echo "All Lean checks passed; generated table is up to date."
