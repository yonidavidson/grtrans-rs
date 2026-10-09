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

pub mod bessel;
pub mod chandra;
pub mod constants;
pub mod four_vector;
pub mod interpolate;
pub mod kerr;
pub mod math;

pub use four_vector::{FourVector, Metric, MINKOWSKI};
