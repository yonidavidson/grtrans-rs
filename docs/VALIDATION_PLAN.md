# Validation plan

The question this plan answers: *how do we know the Rust port computes what
upstream GRTRANS computes, and how do we keep knowing it?*

## 1. Principles

1. Every claim of equivalence is backed by a checked-in fixture and a test.
2. Tolerances are the upstream project's own
   (`run_grtrans_test_problems_public.py`: 1e-2 relative for images and
   spectra, 5% for integrator cross-comparisons, 2e-2 for the toroidal case).
3. Intermediate quantities are validated separately from end-to-end results,
   so failures localize to a module.
4. Where the Rust implementation deliberately differs (LSODA), the deviation
   is measured on the regression problems and the measurement is part of CI.

## 2. Test layers

### Layer 1 — kernel fixtures (bit-level)

Small Fortran drivers (`reference/fortran/*.f90`) call upstream library
routines and dump inputs/outputs to text/raw files. Rust unit tests load the
fixtures and compare:

| Kernel | Comparison | Tolerance |
| --- | --- | --- |
| `zbrent`, `tsum`, `get_weight`/`hunt`/`locate` | exact vs Fortran | 1e-15 rel (identical algorithm) |
| BL/KS metric, LNRF frames, `calc_nullp`, `calcg`, `comoving_ortho` | array comparisons | ~1e-12 rel or documented operation-order tolerance |
| Bessel/Juttner helpers | grid comparison | 1e-10 |
| Chandra table interpolation | grid comparison | 1e-12 |
| polsynch coefficients j,α,ρ (11 outputs) over (ν,B,θ_e,γ_min,p) grid | relative + absolute | 1e-10 rel where smooth; documented near cancellation points |
| `imatrix_4`, `calc_O`, delo single steps | exact algorithms | 1e-13 |

### Layer 2 — ray fixtures (intermediate)

A Fortran debug driver runs single pixels with `debug=1`/`extra=1` and saves
the full per-ray arrays (`x`, `k`, `lambda`, `tpm/tpr`, `rshift`, `s2xi`,
`c2xi`, `j`, `K`, `tau`, `I`). Rust tests compare stage by stage. Covered
pixels:

- THINDISK: center, near photon ring, mid-disk, shadow edge (standard=2).
- FFJET: on-jet, counter-jet, spine/sheath (standard=1, radial).
- HARM: several pixels spanning the flow.
- SPHACC: several impact parameters.

### Layer 3 — end-to-end fixtures

`tests/reference/` runs the full Rust pipeline on the exact inputs in
`reference/fixtures/<case>/inputs.in` and compares to `ivals.f64.bin`:

| Case | Metric | Tolerance |
| --- | --- | --- |
| THINDISK | Σ|ΔI|/Σ|I| over all pixels/Stokes/ν | 1e-2 |
| FFJET | Σ|ΔI|/Σ|I| | 1e-2 |
| HARM | Σ|ΔI|/Σ|I| | 1e-2 |
| SPHACC | profile and spectrum | 1e-1 (upstream uses 10·tol) |
| POWERLAW | Σ|ΔI|/Σ|I| | 2e-2 |
| FFJET delo vs lsoda; FFJET formal vs lsoda | max relative difference of spectra | 5% (upstream) |

Additional internal consistency tests (no Fortran required):

- `λ`-forward vs `λ`-backward integration consistency for delo/formal.
- Invariance checks: $I_\nu/\nu^3$ for thermal emission between frames.
- Stokes-frame rotations leave $I^2-Q^2-U^2-V^2$ and total intensity
  invariant (rotate_emis sanity).
- Flat-space limit (a=0, large r): geodesics reduce to straight lines;
  metric reduces to Minkowski.

## 3. The LSODA deviation

Upstream integrates the transfer equation with ODEPACK LSODA (rtol=1e-6,
atol=1e-8, `hmax=0.1` in affine units for the `[i1:i2]` window). The Rust
port uses an in-house adaptive solver implementing the LSODA algorithm
class (automatic Adams/BDF switching, same tolerances and step limits).

Validation obligations before claiming equivalence:

1. **ODE solver unit tests** against ODEPACK on the actual RHS class
   (piecewise-linear coefficient interpolation): a Fortran driver records
   LSODA solutions for representative (j,K,λ) tables; the Rust solver must
   agree to ≲1e-8 relative at the output points.
2. **Integrator cross-check** on FFJET: delo vs lsoda vs formal spectra must
   agree within the upstream 5% criterion, and the Rust lsoda result must
   match the regenerated Fortran lsoda image within 1e-2.
3. The deviation and its provenance are recorded in `docs/PORTING_MATRIX.md`
   and the crate docs. If test (1) cannot reach the required agreement, the
   fallback is to port the ODEPACK routine itself (tracked as a blocker).

## 4. Fixture regeneration

Fixtures are regenerated only by:

```bash
scripts/build_upstream.sh /path/to/workdir     # fetch + patch + build
GRTRANS_UPSTREAM=/path/to/workdir/upstream \
  python3 scripts/gen_reference.py --outdir reference/fixtures --overwrite
python3 scripts/verify_fixtures.py             # checksums
```

`scripts/verify_fixtures.py` fails if any fixture is missing or its checksum
mismatches its manifest, so tests cannot silently skip.

## 5. CI

`.github/workflows/ci.yml` runs on Linux and macOS:

1. `cargo fmt --check`, `cargo clippy -- -D warnings`.
2. `cargo test` — all layers 1–3 (fixtures are checked in; no Fortran needed
   for the Rust test run).
3. A weekly/`workflow_dispatch` job rebuilds upstream with gfortran and
   regenerates fixtures to detect upstream drift.

## 6. What would falsify equivalence

- A fixture mismatch beyond tolerance in any regression case.
- A kernel fixture mismatch (bit-level layers) that is not explained by an
  operation-order transcription bug.
- Any regression case where the Rust result matches the fixture only after
  tuning constants away from their upstream values.
