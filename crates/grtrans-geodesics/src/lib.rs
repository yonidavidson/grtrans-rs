//! Kerr geodesics for grtrans-rs.
//!
//! The core is a direct translation of geokerr (Dexter & Agol 2009), the
//! semi-analytic Kerr geodesic solver distributed with GRTRANS
//! (`geokerr_wrapper.f` @ c76cb11). On top of it sits the ray assembly from
//! `geodesics.f90`: converting camera pixels to geodesics and returning BL
//! coordinates, wave vectors, affine parameters and turning-point parities.
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

pub mod geokerr;

pub mod camera;
pub mod rays;
