//! Physical constants.
//!
//! Direct translation of `phys_constants.f90` (upstream GRTRANS).
//! All values are reproduced verbatim, including upstream's rounding
//! (e.g. `msun = 1.998e33` rather than the IAU value), because they enter
//! scientific results and must match bit-for-bit behavior.

/// Planck constant (erg s).
pub const H: f64 = 6.626e-27;
/// Boltzmann constant (erg / K).
pub const K: f64 = 1.38e-16;
/// Speed of light (cm / s).
pub const C: f64 = 2.99792458e10;
/// Elementary charge (statC, cgs-Gaussian).
pub const E: f64 = 4.8032e-10;
/// Gravitational constant (cgs).
pub const G: f64 = 6.67e-8;
/// Electron mass (g).
pub const M: f64 = 9.10938188e-28;
/// Proton mass (g).
pub const MP: f64 = 1.67262158e-24;
/// Pi (as evaluated by `acos(-1d0)` in Fortran; equal to the correctly
/// rounded double).
pub const PI: f64 = std::f64::consts::PI;
/// Speed of light squared (cm^2 / s^2).
pub const C2: f64 = C * C;
/// Stefan-Boltzmann constant (cgs).
pub const SIGB: f64 = 5.6704e-5;
/// Solar mass (g); upstream value.
pub const MSUN: f64 = 1.998e33;
/// Thomson cross section (cm^2).
pub const SIGT: f64 = 6.6523e-25;
