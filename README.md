# grtrans-rs

A Rust port of **GRTRANS**, the general-relativistic, polarized radiative
transfer code by [Jason Dexter](https://github.com/jadexter/grtrans)
([Dexter 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.462..115D),
[Dexter & Agol 2009](https://ui.adsabs.harvard.edu/abs/2009ApJ...696.1616D)).

**Status: four of the five upstream regression problems reproduce the
reference results within (and far inside) the upstream tolerances; the HARM
GRMHD reader remains to be ported.** See `docs/PORTING_MATRIX.md` for the
per-module state and `docs/VALIDATION_PLAN.md` for the evidence.

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
kernel fixtures, per-ray intermediate fixtures, and the upstream regression
problems. Reference fixtures were regenerated on this machine and reproduce
the upstream-shipped regression pickles within the upstream thresholds
(`docs/REFERENCE_RESULTS.md`).

| Case | Metric | Upstream tolerance | grtrans-rs |
| --- | --- | --- | --- |
| THINDISK (4 Stokes, 25 ν) | Σ\|ΔI\|/Σ\|I\| | 1e-2 | 1.2e-7 |
| FFJET (POLSYNCHPL, 4 Stokes) | Σ\|ΔI\|/Σ\|I\| | 1e-2 | 4.2e-4 |
| SPHACC (SYNCHTHAV) | Σ\|ΔI\|/Σ\|I\| | 1e-1 | 6.3e-4 image, 5.2e-4 spectrum |
| POWERLAW (POLSYNCHTH, 4 Stokes) | Σ\|ΔI\|/Σ\|I\| | 2e-2 | 7.4e-4 |
| HARM (POLSYNCHTH) | Σ\|ΔI\|/Σ\|I\| | 1e-2 | pending (reader not yet ported) |

Kernel-level fixtures (Bessel functions, Carlson integrals, ZROOTS, GEOMU,
all four emissivity families) agree to ≤1e-13 or better; the geokerr
production path agrees to ≤1e-11 for direct ray quantities. The port also
reproduces upstream's `radtrans_aux` interpolation quirk deliberately, since
it materially shapes the reference images (see `docs/VALIDATION_PLAN.md`
section 3).

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
