# grtrans-rs architecture

## 1. Goals and constraints

1. **Numerical equivalence first.** Every routine that participates in a
   scientific result is a direct translation of the upstream Fortran, with
   the same operation order where it matters for floating-point results
   (accumulations, root finding, interpolation weights).
2. **Layered, testable crates.** Mathematical kernels are independent of I/O
   and language bindings, so each layer can be validated against Fortran
   fixtures.
3. **Reproducibility.** Inputs, outputs, and fixtures are deterministic and
   versioned; no hidden global state (the Fortran used module globals and
   threadprivate state, mirrored in Rust as explicit per-ray structs).
4. **Documented deviations.** Where a literal port is impractical (ODEPACK
   LSODA), the replacement is documented, measured against the original on
   the regression problems, and gated in tests.

## 2. Crate layout

```text
grtrans-rs/
├── Cargo.toml               # workspace
├── crates/
│   ├── grtrans-core/        # constants, math, four-vectors, Kerr metric,
│   │                        # interpolation, Bessel, Chandra table
│   ├── grtrans-geodesics/   # geokerr port, camera geometry, ray sampling
│   ├── grtrans-physics/     # fluid models, emissivities (polsynch etc.)
│   ├── grtrans-transfer/    # polarized transfer integrators (delo/formal/lsoda)
│   ├── grtrans-io/          # input namelists, FITS and binary output
│   ├── grtrans/             # the CLI binary (driver + threading)
│   └── grtrans-python/      # PyO3 bindings (feature-gated)
├── tests/                   # workspace-level integration + reference tests
│   ├── integration/
│   └── reference/           # fixture-driven end-to-end tests
├── reference/fixtures/      # regenerated upstream outputs (see README)
├── scripts/                 # upstream fetch/build, fixture generation
└── docs/
```

Mapping from upstream modules to crates:

| Upstream | Rust |
| --- | --- |
| `phys_constants.f90`, `math.f90`, `class_four_vector.f90`, `kerr.f90`, `interpolate.f90`, `bessel.f90`, `chandra_tab24.f90`, `locate.f`, `hunt.f`, `polint.f`, `polyvl.f` | `grtrans-core` |
| `geokerr_wrapper.f`, `geodesics.f90`, `camera.f90` | `grtrans-geodesics` |
| `emis.f90`, `polsynchemis.f90`, `calc*max*.f90`, `calcgmin.f90`, `fluid.f90`, `fluid_model_*.f90` | `grtrans-physics` |
| `radtrans_integrate.f90`, `rad_trans.f90`, ODEPACK | `grtrans-transfer` |
| `read_inputs.f90`, `fits.f90`, binary output | `grtrans-io` |
| `grtrans_driver.f90`, `pgrtrans.f90`, `grtrans.f90`, `grtrans_program.f90` | `grtrans` |
| `pgrtrans.f90` (f2py interface), `class_geokerr.f90` | `grtrans-python` |

## 3. Type and units conventions

- All floating point is `f64`. Upstream mixed `real(4)` arrays for stored
  simulation data (`fluid_model_*`); where a `real(4)` value feeds the
  physics, the Rust port stores `f32` at the same stage to keep bit-level
  behavior of interpolation comparable, converting to `f64` at use sites.
- **Geometrical units** $G=c=M=1$ inside the metric/geodesic/fluid layer
  (same as upstream).
- **CGS** at the emissivity/transfer layer: `phys_constants` values are
  copied verbatim, including upstream's rounded constants (e.g.
  `msun=1.998e33`, `sigt=6.6523e-25`), because they enter results.
- Metric storage is the upstream 10-component symmetric convention
  `[g00,g01,g02,g03,g11,g12,g13,g22,g23,g33]` in Boyer–Lindquist covariant
  form when produced by `kerr_metric`.
- Four-vectors carry their metric (as upstream `four_Vector`) so that inner
  products work anywhere in the code; `FourVector { data: [f64;4], metric:
  Metric }`.

## 4. Data flow of one ray (mirrors `grtrans_driver.f90`)

```text
camera pixel index i
   │  initialize_camera_geokerr  →  α, β, q², l, u_f, μ_f, s_u, s_m, tpm, tpr
   ▼
geokerr (per pixel)            →  u(τ), μ(τ), t(τ), φ(τ), λ(τ)
   ▼  geodesics assembly       →  x^μ(λ) (BL), k_μ(λ) from calc_nullp
   ▼  fluid model             →  ρ, p, u^μ, b^μ (code units)
   ▼  comoving_ortho          →  g, s2ξ, c2ξ, angle b–k, cosne
   loop over ν, Ṁ, μ:
     convert fluid → cgs emission state (n_e, B, T_e, distribution params)
     calc_emissivity(ν_emit, …) → j_I,Q,U,V, α_I,Q,U,V, ρ_Q,U,V  (11 cols)
     rotate_emis(s2ξ,c2ξ); invariant scaling j·g², α/g
     e%j /= fac; e%K *= L_BH; optical depth τ(s)
     integrate →  I_ν,obs = fac · L_BH · ∫…
     save pixel
```

## 5. Threading model

Upstream uses OpenMP over pixels with `threadprivate` globals. Rust uses
per-ray owned state and `rayon` parallel iteration over pixels; results are
written into a preallocated image array, so bitwise results are identical
to serial execution.

## 6. Determinism

GRTRANS is deterministic: no Monte-Carlo or random sampling in the
scientific path (`PHATDISK` uses a deterministic log-normal distribution
construction). There are therefore no RNG seeds to control. Any future
stochastic feature must use a seeded, documented generator.

## 7. Error handling

Inexact upstream behaviors are preserved rather than "fixed" silently:

- `calc_nullp` clamps negative radicands to zero (upstream `merge`), and
  non-invertible Delo matrices produce `inf`/`NaN` exactly as upstream.
  These cases are asserted in tests.
- Ray failures (geokerr returning no valid points) set the pixel to zero,
  as upstream.
- Rust `Result` is used for I/O and configuration errors only; physics
  kernels are total functions returning floats, matching Fortran.

## 8. Feature flags

| Feature | Effect |
| --- | --- |
| `python` | builds `grtrans-python` PyO3 module (needs Python dev headers) |
| `fits` | FITS output support in `grtrans-io` (pure-Rust `fitsio` crate or hand-rolled writer; see PORTING_MATRIX) |
| `serde` | JSON manifests for fixtures and CLI configs |

## 9. Validation architecture

Three layers, all in CI:

1. **Kernel fixtures** (`tests/unit` + `reference/fortran/`): small Fortran
   drivers dump values of individual routines (metric, calc_nullp, get_weight,
   polsynch coefficients, delo steps); Rust tests compare bit-for-bit or to
   tight tolerances.
2. **Ray fixtures**: full per-ray arrays (x, k, g, j, K, I) for selected
   pixels, produced by a Fortran debug driver; Rust compares intermediate
   stages.
3. **End-to-end fixtures** (`tests/reference`): the six upstream regression
   problems, compared against regenerated upstream outputs with the same
   tolerances as `run_grtrans_test_problems_public.py` (1e-2 relative for
   images, spectra, etc.).
