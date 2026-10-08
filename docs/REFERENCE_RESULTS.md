# Reference results: upstream GRTRANS build and fixture regeneration

This document records the reproducible reference environment, the exact steps
to rebuild upstream GRTRANS, and the results of regenerating the upstream
regression fixtures. Everything here was executed on the machine that hosts
this project on 2026-10-08.

## 1. Reference environment

| Item | Value |
| --- | --- |
| Host | macOS arm64 (Apple Silicon), Darwin 25.x |
| Fortran compiler | `GNU Fortran (Homebrew GCC 15.2.0) 15.2.0` |
| Fortran flags | `-std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -fopenmp -O3 -fPIC` |
| Link flags | `-fopenmp`, cfitsio from Homebrew |
| cfitsio | Homebrew `cfitsio` (libraries at `/opt/homebrew/lib`) |
| Python | CPython 3.14.5 (Homebrew), virtualenv at `~/.venvs/grtrans` |
| NumPy | 2.3.4 |
| Astropy | 8.0.1 |
| Upstream revision | `c76cb11fa1396516f38ba6f972c68fbb5b5ae984` |

## 2. Build recipe

```bash
git clone https://github.com/jadexter/grtrans.git upstream
cd upstream
git checkout c76cb11fa1396516f38ba6f972c68fbb5b5ae984

# 1. apply the documented compatibility patch (namelist declaration order):
python3 /path/to/grtrans-rs/reference/upstream/patch_namelists.py

# 2. configure the build
cat > Makefile.top <<EOF
USEINTEL=0
USEGNU=1
DEBUG=0
XEONPHI=0
PROFILE=0
CFITSIODIR=/opt/homebrew/lib
GRTRANSDIR=$PWD
FFLAGS=-std=legacy -fallow-argument-mismatch
EOF

# 3. build the binary and the static library
make grtrans -j8
make libgrtrans -j8

# 4. build the f2py Python module (used by grtrans_batch.py).
#    The numpy meson backend cannot link the static archive, so build the
#    extension manually:
python -m numpy.f2py -m pgrtrans pgrtrans.f90
# compile pgrtransmodule.c + fortranobject.c (numpy/f2py/src) + wrappers,
# then link with -Wl,-force_load,libgrtrans.a -Wl,-undefined,dynamic_lookup
# (see scripts/build_upstream.sh for the exact commands).

# 5. generate fixtures
GRTRANS_UPSTREAM=$PWD OMP_NUM_THREADS=8 python3 \
    /path/to/grtrans-rs/scripts/gen_reference.py \
    --outdir /path/to/grtrans-rs/reference/fixtures
```

The compatibility patch changes **no semantics**: it only moves `namelist`
statements below the declarations of their members, which `gfortran >= 10`
requires. The patch is stored in `reference/upstream/patch_namelists.py` and
the resulting diff in `reference/upstream/gfortran15.patch`.

## 3. Regenerated fixtures and errors against the shipped pickles

The upstream repository ships regression pickles produced by earlier GRTRANS
runs. Re-running the built binary on the same inputs should reproduce them;
the table below reports the same relative error metric as
`run_grtrans_test_problems_public.py`:
`rel = Σ|ours − pickle| / Σ|pickle|`.

| Case | Fixture (`reference/fixtures/`) | Metric | Relative error | Upstream pass threshold |
| --- | --- | --- | --- | --- |
| THINDISK | `thindisk/ivals` (100×100 px, 4 Stokes, 25 ν) | image | 6.64e-4 | < 1e-2 |
| FFJET (`POLSYNCHPL`) | `ffjet/ivals` (100×100, 4 Stokes, 1 ν) | image | 5.38e-3 | < 1e-2 |
| FFJET unpol (`SYNCHPL`) | `ffjet_unpol/ivals` | image | (no pickle; used for cross-checks) | — |
| FFJET `delo` | `ffjet_delo/ivals`, `spec` | — | compared against lsoda run | Δ < 5% upstream |
| FFJET `formal` | `ffjet_formal/ivals`, `spec` | — | compared against lsoda run | Δ < 5% upstream |
| HARM | `harm/ivals` (150×150, Stokes I) | image | 3.77e-3 | < 1e-2 |
| SPHACC | `sphacc/ivals[:,0,14]` | intensity profile | 2.97e-3 | < 1e-1 (10·tol) |
| SPHACC | `sphacc/spec` | spectrum | 4.54e-3 | < 1e-1 (10·tol) |
| POWERLAW (toroidal) | `powerlaw/ivals` (200×200, 4 Stokes) | image | see manifest | < 2e-2 |

(Exact values are in each fixture's `manifest.json`;
`reference/fixtures/generation_manifest.json` aggregates them.)

Every case regenerates the shipped upstream pickle within the upstream
tolerance, which establishes that (a) the reference environment is faithful
and (b) the shipped pickles are consistent with the recorded revision.

## 4. Fixture provenance and integrity

Every fixture directory contains:

- `inputs.in`, `files.in` — the exact input decks used;
- `*.npy` — NumPy arrays as read from the GRTRANS FITS output;
- `*.f64.bin` — the same arrays as raw little-endian float64, for Rust;
- `manifest.json` — upstream commit, dirty-tree flag, compiler/Python
  versions, input parameters, wall time, SHA-256 checksums of every file,
  and the relative error against the shipped pickle.

`scripts/verify_fixtures.py` re-checks every checksum offline.

## 5. Notes and caveats

- The upstream FITS writer stores float32 (`real(4)`) image data; fixtures
  preserve that, so agreement below ~1e-6 is not expected.
- FFJET `delo`/`formal` runs use 1600 points per ray instead of 400 (as in
  the upstream test) and are cross-compared between integrators; their
  `manifest.json` records inputs and checksums.
- Wall-clock times in manifests are from this machine and are informational
  only.
- The POWERLAW case is the most expensive (200×200×1600 step geodesics);
  if it is absent, `scripts/gen_reference.py --only powerlaw` regenerates it.
