#!/usr/bin/env bash
# Build and run the Fortran kernel-fixture drivers, writing outputs under
# reference/fixtures/fortran/.
#
# Usage: scripts/build_fortran_fixtures.sh WORKDIR
#   WORKDIR must contain a built upstream tree (see scripts/build_upstream.sh).
set -euo pipefail

WORKDIR="${1:?usage: build_fortran_fixtures.sh WORKDIR}"
UPSTREAM="$WORKDIR/upstream"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO_ROOT/reference/fixtures/fortran"
SRC="$REPO_ROOT/reference/fortran"

if [ ! -f "$UPSTREAM/libgrtrans.a" ]; then
    echo "no libgrtrans.a in $UPSTREAM; run build_upstream.sh first" >&2
    exit 1
fi

mkdir -p "$OUT"
cd "$UPSTREAM"

# Build a shared library from the archive (the archive contains a `main`,
# so -force_load into executables would produce duplicate symbols).
if [ ! -f libgrtrans.dylib ]; then
    gfortran -shared -o libgrtrans.dylib \
        -Wl,-force_load,libgrtrans.a \
        -Wl,-install_name,"$PWD/libgrtrans.dylib" \
        -Wl,-undefined,dynamic_lookup \
        -fopenmp -L/opt/homebrew/lib -lcfitsio
fi

for f in "$SRC"/*.f90; do
    name="$(basename "$f" .f90)"
    echo ">> building $name"
    gfortran -std=legacy -fallow-argument-mismatch -O2 -I. \
        "$f" -L. -lgrtrans \
        -L/opt/homebrew/lib -lcfitsio -fopenmp -o "/tmp/$name"
    "/tmp/$name" > "$OUT/$name.txt"
    echo "   wrote $OUT/$name.txt ($(wc -l < "$OUT/$name.txt") lines)"
done
