# Upstream GRTRANS: provenance, licensing, and inventory

This document records the upstream reference implementation used by grtrans-rs,
its licensing, its scientific content, and the exact environment used to build
the reference binary.

## 1. Repository and revision

| Item | Value |
| --- | --- |
| Repository | https://github.com/jadexter/grtrans |
| Clone date | 2026-10-08 |
| Reference revision | `c76cb11fa1396516f38ba6f972c68fbb5b5ae984` |
| Revision date | 2021-03-08 23:06:47 +0100 |
| Commit subject | `testing new intel compiler fix` |
| Branch | `master` |
| Git tags | none |
| Author / maintainer | Jason Dexter (jdexter@mpe.mpg.de) |

All references to "upstream" in this repository mean the revision above.
Line counts in `SOURCE_INVENTORY.md` were generated from that exact revision.

## 2. Authors and scientific references

The code is described in, and should be cited as:

- **Dexter, J. (2016)**, *A public code for general relativistic, polarised
  radiative transfer around spinning black holes*, MNRAS 462, 115
  (the GRTRANS code paper).
- **Dexter, J. & Agol, E. (2009)**, *A fast new public code for computing
  photon orbits in a Kerr spacetime*, ApJ 696, 1616
  (the geokerr algorithm used for geodesics).

Additional physics implemented in the code (from the source and README):

- Shakura & Sunyaev (1973); Novikov & Thorne (1973); Page & Thorne (1974)
  thin disk model (`THINDISK`).
- Dexter & Agol (2011) inhomogeneous "no-zone" disk (`PHATDISK`).
- Bondi (1952) / Michel (1972) spherical accretion (`SPHACC`).
- Broderick & Loeb (2009) semi-analytic jet (`FFJET`).
- Gammie et al. (2003); Noble et al. (2006) HARM GRMHD data (`HARM`).
- Schnittman & Bertschinger (2004); Broderick & Loeb (2006) orbiting hot
  spots (`HOTSPOT`, `SCHNITTMAN`).
- Melrose (1983); Huang et al. (2009); Dexter (2011) polarized synchrotron
  emission/absorption (`POLSYNCHTH`, `POLSYNCHPL`, ...).
- Del Zanna & Bucciantini (2002) integration scheme (`delo`);
  Degl'Innocenti (1985) formal solution (`formal`).
- ODEPACK/LSODA (Hindmarsh; Petzold) for the `lsoda` integrator.
- Walker–Penrose constants and Chandrasekhar (1983) polarization transport;
  Connors, Piran & Stark (1980) ("CPS80") thin-disk polarization degree.
- Beckwith et al. (2008); Shcherbakov & Huang (2011) comoving orthonormal
  frame construction.

## 3. License and attribution

Upstream is distributed under the **MIT License** (file `LICENSE` in the
upstream root):

```
MIT License

Copyright (c) 2016 jadexter

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
```

**Licensing gate: PASSED.** Redistribution and derivative works (including
translation to another language) are permitted provided the copyright notice
and permission notice are preserved. grtrans-rs therefore:

- ships the upstream MIT notice in `LICENSE` alongside the project license,
- preserves attribution in source-file headers where code is a direct
  translation,
- cites Dexter (2016) and Dexter & Agol (2009) in `CITATION.cff` and README.

Third-party code embedded in upstream (all under compatible terms, retained
for provenance):

| Component | Origin | Notes |
| --- | --- | --- |
| `opkda1.f`, `opkda2.f`, `opkdmain.f`, `odepack.f90`, `odepack_aux.f` | ODEPACK (LLNL), LSODA by L. Petzold / A. Hindmarsh | public-domain-style scientific software distributed by LLNL; translated code is subject to validation per `VALIDATION_PLAN.md` |
| `geokerr_wrapper.f`, `polint.f`, `polyvl.f`, `locate.f`, `hunt.f`, `bessel.f90` | Dexter & Agol (2009) geokerr distribution; Numerical Recipes-style utilities | MIT distribution by upstream |
| `bessel.f90` | Incomplete Bessel function routines (NR-style) | MIT distribution by upstream |

No GPL or otherwise incompatible code was identified in the upstream tree.

## 4. Build environment and dependency inventory

Reference build performed for this project (see `docs/REFERENCE_RESULTS.md`
for the full recipe):

| Item | Value |
| --- | --- |
| Platform | macOS (arm64), Darwin 25.x |
| Compiler | GNU Fortran (Homebrew GCC) 15.2.0 |
| Flags | `-std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -fopenmp -O3 -fPIC` |
| External library | cfitsio (Homebrew) for FITS output |
| Python driver | CPython 3.14.5, NumPy 2.3.4, Astropy 8.0.1 |
| Python bindings | `f2py` (manual compile recipe in `docs/REFERENCE_RESULTS.md`) |

Runtime dependencies of upstream:

- **cfitsio** — mandatory for FITS output (`fits.f90`); plain-binary output
  (`cflag=0`) avoids it at runtime but upstream still links it.
- **OpenMP** — optional parallelization of the image loop (`grtrans.f90`,
  `grtrans_driver.f90`); single-thread results are the reference.
- **Python** (optional) — `grtrans_batch.py` driver, `run_grtrans_*.py` tests,
  `pgrtrans.f90` via f2py.
- **Fortran 77/90** — no BLAS/LAPACK, no MPI.

Upstream build patches required with modern compilers (documented, minimal,
and semantically neutral; see `reference/upstream/`):

1. `read_inputs.f90` and 15 fluid-model/program files: move `namelist`
   statements after the type declarations (`gfortran >= 10` requires namelist
   members to be declared before the namelist statement).
2. `Makefile.top` with the local cfitsio path (build configuration only).

No source semantics were changed.

## 5. Scientific functionality inventory

### 5.1 Geodesics

- **geokerr** (Dexter & Agol 2009): semi-analytic Kerr geodesics using
  elliptic integrals/Jacobi elliptic functions. Supports radial tracing
  (`standard=1`, extended emission) and polar-angle tracing to the equatorial
  plane (`standard=2`, thin disks). Camera shapes: `nrotype=1` (polar
  annulus), `=2` (Cartesian image plane).
- Wave four-vector reconstruction from the constants of motion
  (`calc_nullp`), Kerr–Schild and Boyer–Lindquist coordinates, LNRF frames.

### 5.2 Fluid models (upstream names)

| Model | Type | Data source |
| --- | --- | --- |
| `THINDISK` | analytic | NT73/Page–Thorne thin disk |
| `PHATDISK` | analytic-model + precomputed table | Dexter & Agol (2011) |
| `NUMDISK` | numerical data | user-supplied $T_\mathrm{eff}(r,\phi)$ table |
| `HOTSPOT`, `SCHNITTMAN` | analytic | orbiting hot spot |
| `SPHACC` | analytic | Bondi/Michel spherical accretion |
| `FFJET` | numerical data | Broderick & Loeb (2009) dump |
| `HARM`, `HARM3D`, `HARMPI`, `IHARM`, `KORAL`, `KORAL3D`, `THICKDISK`, `MB09` | numerical data | GRMHD dumps (binary) |
| `POWERLAW`, `SARIAF` | analytic | power-law/self-similar models |
| `TOY` | analytic | test model |

### 5.3 Emissivities

`lambda`, `INTERP`, `INTERPPOL`, `BB`, `BBPOL`, `FBB`, `FBBPOL`, `RHO`,
`SYNCHPL`, `SYNCHTH`, `SYNCHTHAV`, `POLSYNCHPL`, `POLSYNCHTH`,
`POLSYNCHSYMTH`, `HYBRIDTH`, `HYBRIDPL`, `HYBRIDTHPL`, `MAXJUTT`, `MAXCOMP`,
`BINS`, `HYBRIDTHBINS`, `SYNCHTHBREMS`, `SYNCHTHAVNOABS`, `BREMS`.
Polarized synchrotron coefficients follow the 11-output convention
`(jI,jQ,jU,jV,αI,αQ,αU,αV,ρQ,ρU,ρV)` (Melrose 1983; Dexter 2011).

### 5.4 Radiative transfer

- `lsoda` — ODEPACK/LSODA adaptive solver (default; linear Stokes 4-vector or
  intensity-only).
- `lsodasph` — experimental spherical Stokes formulation.
- `delo` — Del Zanna & Bucciantini (2002) scheme (4-vector; intensity-only
  falls back to quadrature).
- `formal` — formal solution using the emissivity/opacity operators
  (Degl'Innocenti 1985 style, `calc_O`).
- `quadrature` — intensity-only $\int j e^{-\tau}\,ds$.

### 5.5 Execution modes

- Image mode: camera grid `nn(1)×nn(2)` pixels × `nn(3)` geodesic points;
  loops over $\nu$ (`nfreq`), $\dot M$ (`nmdot`), and observer $\mu$ (`nmu`).
- Planck-convolved model energy spectra (`phatdisk`, `numdisk` tables).
- `extra`/`debug` diagnostics (optical depths, emissivity-weighted averages,
  polarization angles, geodesic debug dumps).
- Output: FITS (cfitsio) or plain binary.

## 6. Known limitations of upstream (as documented in-code)

1. `README`: `EMISTABLE` emissivity is listed under "TO ADD" and is not
   implemented.
2. `lsodasph` spherical-Stokes integration is marked *"in development"* in
   `radtrans_integrate.f90`.
3. `fluid_model_harmpi`/`koral3d`/`iharm`/`mb09`/`thickdisk` have hard-coded
   paths and file-format assumptions for private simulation data; only
   `HARM` (`dump040`), `FFJET` (`m87bl09rfp10xi5a998fluidvars.bin`) data are
   shipped in the repository.
4. The `analytic_pol_rad_trans.py`/tests rely on Python 2/3 compatibility
   shims; some test drivers require the f2py module.
5. No continuous-integration or unit-test harness is present; validation is
   via `run_grtrans_test_problems_public.py` comparison against shipped
   `.p` pickle fixtures.
6. Upstream mixes F77 (`geokerr_wrapper.f`, ODEPACK) and F90 sources; several
   modules use `real(kind=4)` for large arrays (`fluid_model_*`), which the
   Rust port reproduces where it affects results (documented per module).

## 7. Tests and reference outputs shipped upstream

| Artifact | Description |
| --- | --- |
| `test_grtrans_thindisk.p` | THINDISK image `ivals` (100×100 px, 4 Stokes, 25 ν) |
| `test_grtrans_ffjet.p` | FFJET image (100×100, 4 Stokes, 1 ν) for `POLSYNCHPL` |
| `test_grtrans_harm.p` | HARM image (150×150, Stokes I, 1 ν) with `dump040` |
| `test_grtrans_sphacc_intensity.p` | SPHACC 1-D intensity profile (`ivals[:,0,14]`) |
| `test_grtrans_sphacc_spectrum.p` | SPHACC spectrum (25 frequencies) |
| `test_toroidalfield.p` | POWERLAW toroidal-field image (200×200, 4 Stokes) |
| `unit_tests_public.py`, `unit_tests_integration.py`, `simple_radtrans_integrate_tests.py` | auxiliary unit checks (less formal) |
| `run_grtrans_test_problems_public.py` | the canonical regression driver |

grtrans-rs regenerates all of these with its own reference scripts and stores
Rust-friendly fixtures under `reference/fixtures/` (see
`docs/REFERENCE_RESULTS.md`).
