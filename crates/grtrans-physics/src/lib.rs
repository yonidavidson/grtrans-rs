//! Fluid models and emissivities for grtrans-rs.
//!
//! Direct translations of the upstream GRTRANS modules:
//! * `fluid.f90` + `fluid_model_*.f90` -> [`fluid`], [`models`]
//! * `emis.f90` + `polsynchemis.f90`   -> [`emissivity`], [`polsynch`]
// Lints deliberately allowed: the port preserves upstream Fortran
// semantics and constants exactly.
// * `approx_constant`: upstream uses rounded literals (e.g. the
//   Euler-Mascheroni value -0.57721566 in bessel.f90) that must not be
//   replaced by the standard constants.
// * `neg_cmp_op_on_partial_ord`: `!(x > 0.0)` reproduces Fortran
//   `merge(a, b, x > 0)` semantics (NaN takes the false branch), which
//   `x <= 0.0` would not.
// * `too_many_arguments`: upstream routine signatures are preserved.
#![allow(clippy::approx_constant)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)] // index loops mirror Fortran 1-based indexing
#![allow(clippy::chunks_exact_to_as_chunks)]
#![allow(clippy::type_complexity)]
#![allow(clippy::needless_update)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_swap)]
#![allow(clippy::let_and_return)]
#![allow(clippy::items_after_test_module)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::assign_op_pattern)]
#![allow(unused_assignments)] // Fortran re-assignment patterns are preserved

pub mod emissivity;
pub mod fluid;
pub mod models;
pub mod polsynch;
pub mod polsynch_tables;

pub use emissivity::{
    assign_emis_params, calc_emissivity, invariant_emis, rotate_emis, Emis, EmisParams,
};
pub use fluid::{get_fluid_vars, Fluid, FluidArgs, SourceParams};
