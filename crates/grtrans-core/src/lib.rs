//! Core mathematical and geometric kernels for grtrans-rs.
//!
//! This crate is a direct translation of the corresponding upstream GRTRANS
//! (jadexter/grtrans @ c76cb11) modules:
//!
//! * `phys_constants.f90`   -> [`constants`]
//! * `math.f90`             -> [`math`]
//! * `interpolate.f90` + `hunt.f` + `locate.f` -> [`interpolate`]
//! * `class_four_vector.f90` -> [`four_vector`]
//! * `kerr.f90`             -> [`kerr`]
//! * `bessel.f90`           -> [`bessel`]
//! * `chandra_tab24.f90`    -> [`chandra`]
//!
//! Units follow upstream: geometrized units (G = c = M = 1) for the metric,
//! geodesic and fluid layer; CGS at the emissivity layer with the upstream
//! constant values verbatim.

pub mod bessel;
pub mod chandra;
pub mod constants;
pub mod four_vector;
pub mod interpolate;
pub mod kerr;
pub mod math;

pub use four_vector::{FourVector, Metric, MINKOWSKI};
