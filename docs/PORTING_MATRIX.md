# Porting matrix

Status legend: **done** = ported and fixture-validated; **ported** = code
complete, validation pending; **partial** = subset implemented (details);
**blocked** = documented reason; **n/a** = intentionally not ported.

Update this table as work lands. The authoritative per-file line counts are in
`SOURCE_INVENTORY.md`.

## Core

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `phys_constants.f90` | `core::constants` | pending | verbatim constants |
| `math.f90` `tsum`, `dot_product`, `zbrent` | `core::math` | pending | scalar + array Brent |
| `interpolate.f90` (`interp`, `get_weight`, `locate`) | `core::interpolate` | pending | hunt/locate included |
| `class_four_vector.f90` | `core::four_vector` | pending | metric-carrying four-vector |
| `kerr.f90` | `core::kerr` | pending | metric, tetrads, null vectors, PW transport |
| `bessel.f90`, `polint.f`, `polyvl.f` | `core::bessel` | pending | incomplete Bessel |
| `chandra_tab24.f90` + `ch24_vals.txt` | `core::chandra` | pending | table embedded |

## Geodesics

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `INITIALIZE_CAMERA_GEOKERR` | `geodesics::camera::initialize_camera_geokerr` | pending | both `standard` modes, `nrotype` 1/2 |
| `GEOKERR` + `INDEP_MUF` + `GEOMU` + `GEOPHITIME` | `geodesics::geokerr` | pending | elliptic-integral geodesics |
| `SNCNDN`, `ZROOTS`, `LAGUER`, `CALCIMU*`, `ELLPHITMU` | `geodesics::geokerr` support | pending | |
| `geodesics.f90` ray assembly | `geodesics::rays` | pending | λ, x, k, tpm/tpr |
| `camera.f90` | `geodesics::camera` | pending | FITS header handling delegated to io |

## Emissivity

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `emis.f90` dispatch (+BB, FBB, BBPOL, RHO, INTERP, hybrids) | `physics::emissivity` | pending | |
| `polsynchemis.f90` thermal (`POLSYNCHTH`, `SYNCHTH`, `SYNCHTHAV`) | `physics::polsynch::thermal` | pending | |
| `polsynchemis.f90` power law (`POLSYNCHPL`, `SYNCHPL`) | `physics::polsynch::powerlaw` | pending | |
| `calc_maxjutt.f90`, `calc_maxcomp.f90`, `calcgmin.f90` | `physics::polsynch` helpers | pending | |
| `rotate_emis`, `invariant_emis` | `physics::emissivity` | pending | |

## Fluid models

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `fluid.f90` dispatch/geometry | `physics::fluid` | pending | per-model source params |
| `THINDISK` | `physics::models::thindisk` | pending | |
| `PHATDISK` | `physics::models::phatdisk` | pending | |
| `NUMDISK` | `physics::models::numdisk` | pending | data-file model |
| `HOTSPOT`, `SCHNITTMAN` | `physics::models::hotspot` | pending | |
| `SPHACC` | `physics::models::sphacc` | pending | radial table included in source upstream |
| `FFJET` | `physics::models::ffjet` | pending | uses shipped binary dump |
| `POWERLAW`, `SARIAF`, `TOY` | `physics::models::{powerlaw,sariaf,toy}` | pending | |
| `HARM` | `physics::models::harm` | pending | uses shipped `dump040` |
| `HARM3D`, `HARMPI`, `IHARM`, `KORAL`, `KORAL3D`, `MB09`, `THICKDISK` | `physics::models::*` | pending | private-data formats; may be blocked (documented) if fixtures cannot be produced |

## Transfer

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `radtrans_integrate.f90` delo | `transfer::delo` | pending | incl. thin-step branch |
| `radtrans_integrate.f90` formal | `transfer::formal` | pending | `calc_O` eigenvalues |
| `radtrans_integrate.f90` quadrature | `transfer::quadrature` | pending | nvals=1 fallback |
| LSODA (`opkda*.f`) linear/nonlinear Stokes RHS+Jacobian | `transfer::lsoda` | pending | port or validated substitute; see VALIDATION_PLAN §3 |
| spherical Stokes (`lsodasph`) | `transfer::lsoda_sph` | pending | marked experimental upstream |

## Driver, IO, bindings

| Upstream | Rust item | Status | Notes |
| --- | --- | --- | --- |
| `read_inputs.f90` | `io::inputs` | pending | namelist parser + writer |
| `fits.f90` | `io::fits` | pending | output keys matching upstream |
| plain binary output | `io::binary` | pending | README format |
| `grtrans_driver.f90` | `grtrans::driver` | pending | |
| `grtrans.f90`/`grtrans_program.f90` | `grtrans::main` | pending | `files.in` entry |
| `pgrtrans.f90`, `class_geokerr.f90` | `grtrans-python` | pending | PyO3 |
| Python test harness | `tests/reference` + `scripts/` | pending | fixture regeneration |

## Intentional omissions (n/a)

| Upstream | Reason |
| --- | --- |
| `pgriter.py`, `pgrface.py`, `ppslave.py`, `commit_grtrans.py` | cluster/local workflow scripts, not scientific |
| `fits_prealwin_07102014.f90`, `emis_old.f90` | historical variants, unreferenced by the build |
| `geokerr_interface.py`, development scripts | superseded by the Rust test harness |
