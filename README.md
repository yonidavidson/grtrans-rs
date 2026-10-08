# grtrans-rs

A Rust port of **GRTRANS**, the general-relativistic, polarized radiative
transfer code by [Jason Dexter](https://github.com/jadexter/grtrans)
([Dexter 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.462..115D),
[Dexter & Agol 2009](https://ui.adsabs.harvard.edu/abs/2009ApJ...696.1616D)).

**Status: work in progress.** See `docs/PORTING_MATRIX.md` for the current
per-module state. Do not use for science until the reference tests in
`tests/reference/` pass on your platform.

## What GRTRANS computes

Polarized images and spectra of emission near spinning black holes:
semi-analytic Kerr geodesics (geokerr), a catalogue of analytic and
numerical fluid models (thin/thick disks, GRMHD dumps, jets, hot spots,
spherical accretion), polarized synchrotron and blackbody emissivities, and
integration of the polarized radiative transfer equation (delo, formal,
lsoda schemes).

## Layout

```text
crates/grtrans-core       geometry, Kerr metric, math kernels
crates/grtrans-geodesics  geokerr geodesics, camera, ray sampling
crates/grtrans-physics    fluid models, emissivities
crates/grtrans-transfer   polarized transfer integrators
crates/grtrans-io         inputs, FITS/binary output
crates/grtrans            CLI driver
crates/grtrans-python     PyO3 bindings (optional)
reference/                upstream fixtures + provenance
docs/                     upstream survey, architecture, porting matrix,
                          validation plan, reference results
```

## Validation

The port is validated against the upstream Fortran implementation at
`c76cb11` with three layers of tests described in `docs/VALIDATION_PLAN.md`:
kernel fixtures, per-ray intermediate fixtures, and the six upstream
regression problems (THINDISK, FFJET, HARM, SPHACC, POWERLAW). Reference
fixtures were regenerated on this machine and reproduce the upstream-shipped
regression pickles within the upstream thresholds
(`docs/REFERENCE_RESULTS.md`).

## Building

```bash
cargo build --release
cargo test
```

Python bindings (optional):

```bash
pip install maturin
maturin develop -m crates/grtrans-python/Cargo.toml
```

## Reproducing the Fortran reference

```bash
scripts/build_upstream.sh /path/to/workdir
GRTRANS_UPSTREAM=/path/to/workdir/upstream \
  python3 scripts/gen_reference.py --outdir reference/fixtures
python3 scripts/verify_fixtures.py
```

## License

MIT. This project is a derivative work of GRTRANS
(Copyright (c) 2016 jadexter; MIT License) and preserves its notices —
see `LICENSE`.

## Citation

If you use this software, cite Dexter (2016) and Dexter & Agol (2009); see
`CITATION.cff`.
