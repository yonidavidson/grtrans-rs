//! Kerr geodesics for grtrans-rs.
//!
//! The core is a direct translation of geokerr (Dexter & Agol 2009), the
//! semi-analytic Kerr geodesic solver distributed with GRTRANS
//! (`geokerr_wrapper.f` @ c76cb11). On top of it sits the ray assembly from
//! `geodesics.f90`: converting camera pixels to geodesics and returning BL
//! coordinates, wave vectors, affine parameters and turning-point parities.

pub mod geokerr;

pub mod camera;
pub mod rays;
