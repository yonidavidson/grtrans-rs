//! Fluid models and emissivities for grtrans-rs.
//!
//! Direct translations of the upstream GRTRANS modules:
//! * `fluid.f90` + `fluid_model_*.f90` -> [`fluid`], [`models`]
//! * `emis.f90` + `polsynchemis.f90`   -> [`emissivity`], [`polsynch`]

pub mod emissivity;
pub mod fluid;
pub mod models;
pub mod polsynch;

pub use emissivity::{
    assign_emis_params, calc_emissivity, invariant_emis, rotate_emis, Emis, EmisParams,
};
pub use fluid::{get_fluid_vars, Fluid, FluidArgs, SourceParams};
