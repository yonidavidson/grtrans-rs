#!/usr/bin/env bash
# Build the upstream Fortran GRTRANS reference and its f2py Python module.
#
# Usage: scripts/build_upstream.sh WORKDIR [CFITSIO_PREFIX]
#   WORKDIR        directory that will contain ./upstream (created if needed)
#   CFITSIO_PREFIX prefix containing lib/libcfitsio (default: /opt/homebrew/lib)
#
# Result: WORKDIR/upstream/grtrans (binary), libgrtrans.a, pgrtrans*.so
#
# Tested with: gfortran 15.2.0 (Homebrew), cfitsio (Homebrew), CPython 3.14,
# numpy 2.3.4 on macOS arm64.  Requires: git, make, gfortran, python3 with
# numpy (and meson/ninja if the f2py module is built through numpy's backend).
set -euo pipefail

WORKDIR="${1:?usage: build_upstream.sh WORKDIR [CFITSIO_PREFIX]}"
CFITSIO_PREFIX="${2:-/opt/homebrew/lib}"
UPSTREAM_COMMIT="c76cb11fa1396516f38ba6f972c68fbb5b5ae984"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

if [ ! -d upstream ]; then
    git clone https://github.com/jadexter/grtrans.git upstream
fi
cd upstream
git fetch origin || true
git checkout "$UPSTREAM_COMMIT"

echo ">> applying gfortran>=10 namelist compatibility patch"
python3 "$REPO_ROOT/reference/upstream/patch_namelists.py"

cat > Makefile.top <<EOF
USEINTEL=0
USEGNU=1
DEBUG=0
XEONPHI=0
PROFILE=0
CFITSIODIR=$CFITSIO_PREFIX
GRTRANSDIR=$PWD
FFLAGS=-std=legacy -fallow-argument-mismatch
EOF

echo ">> building grtrans binary and library"
make grtrans -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
make libgrtrans -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"

# --- f2py module -----------------------------------------------------------
# numpy's meson backend cannot link a static Fortran archive properly on
# macOS, so we generate the wrappers with f2py and compile/link manually.
PYTHON="${GRTRANS_PYTHON:-python3}"
if ! "$PYTHON" -c "import numpy" 2>/dev/null; then
    echo ">> python numpy not available; skipping pgrtrans module" >&2
    exit 0
fi

echo ">> building pgrtrans f2py module"
"$PYTHON" -m numpy.f2py -m pgrtrans pgrtrans.f90

PYINC="$("$PYTHON" -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
NPINC="$("$PYTHON" -c 'import numpy; print(numpy.get_include())')"
F2SRC="$("$PYTHON" -c 'import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))')"
EXTSUFFIX="$("$PYTHON" -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"

gfortran -c -O2 -fPIC -std=legacy -fallow-argument-mismatch -I. \
    pgrtrans-f2pywrappers2.f90 -o pgrtrans-f2pywrappers2.o
cc -c -O2 -fPIC -I"$PYINC" -I"$NPINC" -I"$F2SRC" \
    pgrtransmodule.c -o pgrtransmodule.o
cc -c -O2 -fPIC -I"$PYINC" -I"$NPINC" -I"$F2SRC" \
    "$F2SRC/fortranobject.c" -o fortranobject.o
gfortran -shared -o "pgrtrans$EXTSUFFIX" \
    pgrtransmodule.o fortranobject.o pgrtrans-f2pywrappers2.o \
    -Wl,-force_load,libgrtrans.a -Wl,-undefined,dynamic_lookup \
    -L"$CFITSIO_PREFIX" -lcfitsio -fopenmp

echo ">> done. binary: $PWD/grtrans  module: $PWD/pgrtrans$EXTSUFFIX"
