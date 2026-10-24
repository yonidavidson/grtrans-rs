# Porting matrix

Status legend: **done** = ported and fixture-validated; **ported** = code
complete, validation pending; **partial** = subset implemented (details);
**blocked** = documented reason; **n/a** = intentionally not ported.

Per-file line counts: `docs/SOURCE_INVENTORY.md`. Validation evidence:
`docs/VALIDATION_PLAN.md` and `docs/REFERENCE_RESULTS.md`.

## Core (all validated against Fortran fixtures)

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `phys_constants.f90` | `core::constants` | done | verbatim constants |
| `math.f90` `tsum`, `zbrent` | `core::math` | done | scalar + array Brent |
| `interpolate.f90`, `hunt.f`, `locate.f` | `core::interpolate` | done | f64 + f32 variants (overload semantics) |
| `class_four_vector.f90` | `core::four_vector` | done | metric-carrying four-vector |
| `kerr.f90` | `core::kerr` | done | metric (f64/f32), frames, `calc_nullp`, comoving orthonormal frame, Walker–Penrose transport, polarization vectors, plunging/rms velocities, surface integrals |
| `bessel.f90`, `polint.f`, `polyvl.f` | `core::bessel` | done | bit-level fixture agreement (1e-14) |
| `chandra_tab24.f90` | `core::chandra` | done | f64 and f32 variants; table embedded |

## Geodesics (validated against the production path)

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `INITIALIZE_CAMERA_GEOKERR` | `geodesics::geokerr::camera` | done | standard 1/2, `nrotype` 1/2 |
| `GEOKERR`, `INDEP_MUF`, `GEOMU`, `GEOR`, `GEOPHITIME` | `geodesics::geokerr` | done | camera 1e-12, direct ray quantities 1e-11; λ/t cancellation-limited (see plan) |
| Carlson integrals, Jacobi functions, ZROOTS/LAGUER, GAULEG | `geodesics::geokerr::{elliptic,special}` | done | RF/RC/RD/RJ 1e-14; quartic roots 1e-13 |
| `geodesics.f90` ray assembly | `geodesics::rays` | done | matches production `initialize_geodesic` |

Documented upstream quirks reproduced or resolved:

- `standard=2` leaves `TPRARR` unassigned; the port uses 0 (validated:
  standard=2 references are insensitive to it).
- `standard=1` leaves `TPMARR`/`MUFARR` unassigned; the port uses 0 (and the
  fixture driver zeroes them for determinism).
- The MUFILL path writes up to KEXT elements past the output arrays
  upstream; the port guards these writes (values are never read).
- `geodesics.f90` multiplies `phi0` by the f32 value of π
  (`acos(-1.)`); reproduced exactly.

## Emissivity

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `emis.f90` dispatch, BB/FBB/BBPOL, lambda | `physics::emissivity` | done | THINDISK validated end-to-end at 1.2e-7 |
| `polsynchemis.f90` `polsynchpl` | `physics::polsynch` | done | matches Fortran to all printed digits on the fixture grid |
| `polsynchemis.f90` `synchpl` | `physics::polsynch` | done | columns 1/5 compared (others uninitialized upstream) |
| `bnu` | `physics::polsynch::bnu` | done | |
| `polsynchemis.f90` `polsynchth` (thermal) | `physics::polsynch` | pending | needed for HARM |
| `rotate_emis`, `invariant_emis` | `physics::emissivity` | done | |
| `INTERP`/`INTERPPOL`, hybrids, MAXJUTT/MAXCOMP, BREMS | — | blocked | table data not shipped; documented |

## Fluid models

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `fluid.f90` dispatch/geometry | `physics::fluid` | done (ported models) | explicit state instead of module globals |
| `THINDISK` | `physics::models::thindisk` | done | end-to-end 1.2e-7 vs reference |
| `FFJET` | `physics::models::ffjet` | done | fluid fixture at f32 precision; end-to-end within the measured LSODA deviation |
| `SPHACC` | `physics::models::sphacc` | done | stored profile tables embedded; end-to-end image 1.8e-2 / spectrum 7.0e-3 |
| `SPHACC` | `physics::models::sphacc` | done | end-to-end image 6.3e-4 / spectrum 5.2e-4 |
| `POWERLAW` | `physics::models::powerlaw` | done | end-to-end 7.4e-4 (nrotype=1 circular camera validated here too) |
| `SARIAF` | `physics::models::sariaf` | done | self-contained; fluid fixture 2.3e-6; end-to-end 4.1e-3 (upstream LSODA error on thick pixels, VALIDATION_PLAN §3.2) |
| `PHATDISK`, `HOTSPOT`, `SCHNITTMAN`, `TOY` | — | pending | self-contained analytic, not yet ported |
| `NUMDISK` | — | blocked | requires external data files not shipped |
| `HARM` | `physics::models::harm` | done | formatted dump reader, MKS→KS→BL transforms, LNRF velocities, 2D interpolation (single snapshot); fluid fixture 7.6e-6; end-to-end 4.8e-2 with the upstream integrator deviation documented in VALIDATION_PLAN §3.1 |
| `HARM3D`, `HARMPI`, `IHARM`, `KORAL*`, `MB09`, `THICKDISK` | — | blocked | require simulation dumps not shipped with the upstream repository |

## Transfer

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `radtrans_integrate.f90` delo | `transfer::radtrans_integrate_delo` | done | exact scheme; agrees with the exact ODE solution to ~0.4% on FFJET |
| `radtrans_integrate.f90` formal | `transfer::radtrans_integrate_formal` | done | `calc_O` matrix exponential |
| quadrature | `transfer::radtrans_integrate_quadrature` | done | |
| LSODA (`opkda*.f`, 28k lines) | `transfer::lsoda` | **substituted** | Four-Stokes: stiff-stable exponential integrator (exact frozen-coefficient update with adaptive step doubling); matches an independent exact ODE solution to 3e-6 on thick pixels (VALIDATION_PLAN §3.2). Intensity-only (`nequations==1`): unconditionally stable exponential integrator with bounded adaptive substepping (the scalar equation is stiff; upstream LSODA switches to BDF). Both reproduce upstream's `radtrans_aux` one-index value-interpolation quirk deliberately (see VALIDATION_PLAN §3); end-to-end agreement is ≤8e-4 on all regression cases. A faithful ODEPACK port remains an open item. |
| spherical Stokes (`lsodasph`) | — | blocked | experimental upstream ("in development") |

## Driver, IO, bindings

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `grtrans_driver.f90` | `grtrans::driver` | done | single-point and integrated paths validated (THINDISK, FFJET ray chain) |
| `read_inputs.f90` | `grtrans::cli` (in progress) | partial | namelist parsing |
| FITS output (`camera.f90` `kwrite_raytrace_camera`, cflag=1) | `grtrans::fits::write_fits` | done | dependency-free writer (2880-byte blocks, big-endian f32); primary + one IMAGE extension per frequency, NX/NY/NU at positional cards 6/7/8 as the upstream Python reader assumes. Validated with astropy + the upstream reader logic on THINDISK (1.14e-7, ab exact) and a structural round-trip test |
| plain binary output | — | pending | |
| `grtrans.f90`/`grtrans_program.f90` | `grtrans` CLI | partial | orchestrator in place, `files.in` entry pending |
| `pgrtrans.f90` | `grtrans-python` | pending | PyO3 skeleton only |

## Intentional omissions (n/a)

| Upstream | Reason |
| --- | --- |
| `pgriter.py`, `pgrface.py`, `ppslave.py`, `commit_grtrans.py` | cluster/local workflow scripts |
| `fits_prealwin_07102014.f90`, `emis_old.f90` | historical variants, unreferenced by the build |
| development Python scripts | superseded by the Rust test harness |
