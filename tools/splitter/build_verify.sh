#!/bin/bash
# Build the C++ splitter translation, run it on the bundled test case, and
# byte-compare its outputs (stdout + the split grid) against the gfortran
# reference goldens in ref/. Self-contained: works from a fresh clone.
# Requires only g++ (C++17).
#
# splitter is the CFL3D block/input-file splitter (NASA Langley v6.7). It reads
# a control deck on stdin plus the unsplit grid (grdflat5.inp / grdflat5.fmt)
# and writes the split, multigridable grid grdflat5.bin (Fortran unformatted).
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
CPP="$ROOT/cpp"; REF="$ROOT/ref"; RUN="$ROOT/run"

echo "=== compile + link ==="
g++ -std=c++17 -O2 -ffp-contract=off -I"$CPP" "$CPP"/*.cpp "$CPP"/runtime/*.cpp -o "$CPP/splitter"
xattr -c "$CPP/splitter" 2>/dev/null || true
codesign -s - -f "$CPP/splitter" 2>/dev/null || true
echo "LINK OK: $(ls -la "$CPP/splitter" | awk '{print $5}') bytes"

echo "=== run (fresh sandbox) ==="
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$REF/grdflat5.inp" "$REF/grdflat5.fmt" "$REF/split.inp_1blk" "$RUN/"
( cd "$RUN" && "$CPP/splitter" < split.inp_1blk > main.txt 2>err.txt ) || echo "(nonzero exit)"

rc=0
echo "=== compare vs golden ==="
if cmp -s "$RUN/main.txt" "$REF/main.txt.golden"; then
  echo "stdout (main.txt)     : IDENTICAL"
else
  echo "stdout (main.txt)     : DIFFERS"; rc=1
fi
if cmp -s "$RUN/grdflat5.bin" "$REF/grdflat5.bin.golden"; then
  echo "split grid (grdflat5.bin): BYTE-IDENTICAL"
else
  echo "split grid (grdflat5.bin): DIFFERS"; rc=1
fi
[ "$rc" = 0 ] && echo "RESULT: BYTE-EXACT MATCH (stdout + split grid)" \
             || echo "RESULT: MISMATCH (see above)"
exit $rc
