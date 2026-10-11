# Progress log and handoff notes

This file records the state of the port, what is validated, and the exact
next steps. It is updated at each milestone.

## Current state (see docs/VALIDATION_PLAN.md for evidence)

Validated end-to-end against regenerated upstream fixtures (Fortran binary
at `jadexter/grtrans` @ `c76cb11`):

| Case | Upstream tolerance | grtrans-rs result |
| --- | --- | --- |
| THINDISK (BBPOL, 4 Stokes, 25 ν) | 1e-2 | 1.2e-7 |
| FFJET (POLSYNCHPL, 4 Stokes) | 1e-2 | 4.2e-4 |
| SPHACC (SYNCHTHAV, image/spectrum) | 1e-1 | 6.3e-4 / 5.2e-4 |
| POWERLAW (POLSYNCHTH, 4 Stokes) | 2e-2 | 7.4e-4 |
| HARM (POLSYNCHTH) | 1e-2 | **not yet ported (next task)** |

Kernel and per-ray fixtures: Bessel (1e-14), Carlson (1e-14), ZROOTS
(1e-13), GEOMU (1e-13), all four emissivity families (≤8.7e-10), geokerr
production path (camera 1e-12, direct ray quantities 1e-11), FFJET fluid
(f32 level), FFJET per-ray geodebug chain.

All 26 cargo test targets pass. The CLI reproduces THINDISK through the full
file→pipeline→binary path (1.1e-7); the PyO3 wheel builds with maturin and
reproduces THINDISK (1.1e-7).

## Key upstream quirks reproduced (do not "fix" them)

1. `radtrans_aux` (radtrans_integrate.f90) interpolates the transfer
   coefficients one sample deeper than the λ bracket. This materially
   shapes the reference images; the port reproduces it in both integrators
   (`transfer::lsoda`) and includes one extra sample beyond the trimmed
   window so the shifted interpolation can reach zero-emission samples.
   Reproducing it improved FFJET 1.7e-2→4.2e-4, SPHACC 1.8e-2→6.3e-4,
   POWERLAW 5.1e-2→7.4e-4.
2. Fortran NaN comparison semantics: `if (max.gt.errtol) goto 1` exits the
   Carlson convergence loops on NaN; the port uses
   `if !(max > errtol) { break; }` in `geokerr::elliptic`. Without this,
   rays with `uf == u2` (a common MUFILL endpoint) loop forever in RF.
3. Unassigned camera arrays: standard=2 leaves TPRARR unassigned; standard=1
   leaves TPMARR/MUFARR unassigned. The port uses 0 deterministically.
4. MUFILL writes up to KEXT elements past the output arrays upstream; the
   port guards those writes (never read upstream).
5. `geodesics.f90` multiplies phi0 by the f32 value of π (`acos(-1.)`);
   reproduced.
6. `get_weight` reads out of bounds for x outside the table; the port clamps
   to the edge interval (f32 variant; affects only innermost SPHACC points).
7. Fortran table literals with a leading decimal point (`.13733`) must be
   parsed as 0.13733 (`scripts/extract_sphacc_tables.py` handles this; a
   naive regex once produced 13733 and a 6500x error).
8. `synchpl`/`synchemis` leave columns 2-4, 6-11 of their 11-column output
   uninitialized; only columns 1/5 are meaningful.

## Next task: HARM

Goal: reproduce `reference/fixtures/harm/ivals.f64.bin` within 1e-2.

Inputs (from `reference/fixtures/harm/manifest.json`; check it):
`fname='HARM'`, `ename='POLSYNCHTH'`, `nvals=1`, standard=1, nn=[150,150,400],
gridvals=[-15,15,-15,15], uout=0.01, mbh=6.2e8, spin=0.9375, mu=0.6428,
nfreq/fmin/fmax per manifest, `hhfile='dump040'`, `hdfile='dump'`,
`findf=40`, `gmin=1`, `muval=0.25` (check all in the manifest).

The 62 MB `dump040` (formatted, 34 columns × 65536 rows) is in
`~/dev/private/grtrans-rs-work/upstream/dump040`; copy it into
`reference/fixtures/harm/` (like the FFJET binary) so tests are
self-contained.

### Pieces to port (in `fluid_model_harm.f90` + `fluid.f90`)

1. `transformbl2mksh(r, th, x1, x2, h)` (lines 58-71) — BL→MKS coordinate
   transform used by `harm_vals`; note `h` (the MKS flattening) comes from
   the header (header(11) is `gam`; check which header element is `h`:
   `read_harm_data_header` sets `h=header(nhead-1)` with nhead=26 →
   header(25)).
2. `read_harm_data_header` (275-298): reads the first line of `hfile`
   (26 values): nx1, nx2, startx1, startx2, dx1, dx2, asim, gam, h.
3. `read_harm_data_file` (300-408): reads `dfile` + 3-digit index
   (`'dump'//'040'` = `dump040`); first line = 26-value header; then
   `(nx1*nx2)/nx2 = nx1` blocks of `nx2` lines × 34 values. Layout per
   line: cols 1-4 = (x1, x2, r, th) [check: `grid=transpose(data(1:4,:))`,
   `x1_arr=grid(:,1)`, `x2_arr=grid(:,2)`, `r_arr=grid(:,3)`,
   `th_arr=grid(:,4)`], col 5 = rho, col 6 = p, cols 14-17 = u^t,u^r,u^θ,u^φ
   (MKS), cols 22-25 = b (MKS), col 34 = gdet.
   Then: `uks=umksh2uks(u, r_arr, x2_arr, h)` (defined in the same file,
   lines 73-98: ur = r*u^r_mks; dthdx2 = π(1+(1-h)cos(2π x2)); uth =
   u^θ_mks*dthdx2), `u = uks2ubl(uks, r_arr, asim)` (already in
   `grtrans-core::kerr`), same for b, then assign the f64 BL metric.
4. `harm_vals` (100-265): 2D interpolation.
   - `uniqx2 = x2_arr(1:nx2)` (contiguous); `uniqx1 = x1_arr(1:(nx2-1)*nx1+1:nx2)`
     (strided by nx2!); `uniqr = exp(uniqx1)`;
     `uniqth = π*uniqx2 + (1-h)/2*sin(2π uniqx2)`.
   - `transformbl2mksh(zr, theta, x1, x2, h)` per point.
   - `lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1`;
     `lx2` similarly; clamps: `lx2=max(lx2,1)`; `ux2=min(lx2+1,nx2)`;
     (note: lx1 is NOT clamped upstream — rays outside the grid read
     out-of-bounds; clamp deterministically in Rust and document).
   - pole handling: `dth = uniqth(ux2)-uniqth(lx2)` if ux2≠lx2 else
     `uniqth(1)`; `td=|theta-uniqth(lx2)|/dth`;
     `rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))`; if
     `uniqr(lx1) <= 1+sqrt(1-a²)` then `rd=1` (nearest neighbour inside
     the horizon-adjacent zone).
   - bilinear interpolation over the 2D indices
     `x1l=(lx1-1)*nx2`, `x1u=(ux1-1)*nx2`, `x2l=lx2-1`, `x2u=ux2-1`.
   - time interpolation: for nt=1 this reduces to the single snapshot
     (check `tindx`/`ttd` logic; `tt = bl2ks(r, t, a) - bl2ks(r(1),0,a)`,
     `tt = -tt + 1e-8`, `tstep` from the data).
5. `get_harm_fluidvars` (fluid.f90:689): calls `harm_vals(x0,a,f%rho,f%p,
   f%b,f%u,f%bmag)`.
6. `convert_fluidvars_harm` (fluid.f90 ~957-973): `scale_sim_units` for
   rho/p/b/u depending on `sp%gminval`; for `gminval>=1` uses `monika_e`
   with `rlow=1/muval-1`, `rhigh=gminval*(1/muval-1)` (same as POWERLAW);
   for gminval<0 the nonthermal b2 model. Read the exact code first.
7. Test: mirror `crates/grtrans/tests/reference_sphacc.rs` (or
   reference_powerlaw.rs for parallelism); fixture layout (npts, nvals,
   nfreq) per the manifest; tolerance 1e-2.

### Pitfalls learned on similar models

- The dump is formatted text; parse with `str::split_whitespace` +
  `f64::parse` (Fortran may write `D` exponents — handle `d`/`D`).
- Preserve f32/f64 narrowing points exactly (declare where upstream's
  variables are `real` vs `real(8)`).
- `harm_vals` arrays `ppi/rhoi/...` are `(size(x0), 2**(ndim+1))` with
  `ndim=2` → 8 corner values; check how they are filled (the loop after the
  index computation) before writing the interpolation.
- Commit after the reader works (a `test_harm_fluid.f90` fixture comparing
  rho/p/b/u at a few hundred points is the fastest way to debug, like
  `test_ffjet_fluid.f90`).

## Then (optional, in order of value)

1. A faithful ODEPACK LSODA port (bit-exactness; currently ≤8e-4 agreement
   via the substitute).
2. FITS output (`camera.f90` `kwrite_raytrace_camera` cflag=1) — cfitsio
   dependency; binary output already implemented.
3. Remaining fluid models requiring external data (PHATDISK, NUMDISK,
   HOTSPOT, SCHNITTMAN, SARIAF, TOY) — blocked without data files.
4. `lsodasph` (experimental upstream, marked "in development").

## How to run everything

```bash
source ~/.cargo/env
cd ~/dev/private/grtrans-rs
cargo test --release --workspace --exclude grtrans-python   # all layers
cargo clippy --workspace --exclude grtrans-python --all-targets  # clean

# Fortran reference regeneration (needs the built upstream):
bash scripts/build_fortran_fixtures.sh ~/dev/private/grtrans-rs-work
python3 scripts/verify_fixtures.py

# Python bindings:
~/.venvs/grtrans/bin/maturin build --release -m crates/grtrans-python/Cargo.toml --out /tmp/wheels
~/.venvs/grtrans/bin/pip install --force-reinstall /tmp/wheels/*.whl
~/.venvs/grtrans/bin/python tests/python/test_python_api.py
```

Upstream work tree: `~/dev/private/grtrans-rs-work/upstream` (patched for
gfortran≥10 namelist placement; `libgrtrans.dylib`, `grtrans` binary, f2py
modules built).
