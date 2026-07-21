#!/bin/bash
# Build the C++ ronnie translation, run it on the test case, and byte-compare
# its outputs against the gfortran reference (goldens in ref/). Self-contained:
# works from a fresh clone. Requires only g++ (C++17).
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
CPP="$ROOT/cpp"; REF="$ROOT/ref"; RUN="$ROOT/run"

echo "=== compile + link ==="
cd "$CPP"
rm -f ./*.o ronnie
g++ -std=c++17 -O2 -ffp-contract=off *.cpp -o ronnie
# macOS: clear quarantine + ad-hoc sign (no-op elsewhere)
xattr -c ronnie 2>/dev/null || true
codesign -s - -f ronnie 2>/dev/null || true
echo "LINK OK: $(ls -la ronnie | awk '{print $5}') bytes"

echo "=== run (fresh sandbox) ==="
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$REF/grid.bin" "$REF/ronnie.inp" "$RUN/"
cd "$RUN"
"$CPP/ronnie" < ronnie.inp > run.stdout 2>run.stderr || echo "(nonzero exit)"

rc=0
echo "=== compare vs golden ==="
if cmp -s patch.bin  "$REF/patch.bin.golden";  then echo "patch.bin : IDENTICAL"; else echo "patch.bin : DIFFERS"; rc=1; fi
if cmp -s ronnie.out "$REF/ronnie.out.golden"; then echo "ronnie.out: IDENTICAL"; else echo "ronnie.out: DIFFERS"; rc=1; fi
[ "$rc" = 0 ] && echo "RESULT: BYTE-EXACT MATCH" || echo "RESULT: MISMATCH (see above)"
exit $rc
