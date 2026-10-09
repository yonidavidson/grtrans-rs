//! geokerr: semi-analytic Kerr geodesics by Dexter & Agol (2009).
//!
//! Direct translation of `geokerr_wrapper.f` (upstream GRTRANS @ c76cb11),
//! which is self-contained Fortran 77. The module is split into:
//!
//! * [`elliptic`] - Carlson elliptic integrals (RF, RC, RD, RJ) and the
//!   cubic/quartic integral functions for the radial motion.
//! * [`special`]  - Gauss-Legendre nodes, polynomial root finding (ZROOTS,
//!   LAGUER) and Jacobi elliptic functions (SNCNDN).
//! * [`imu`]      - mu-motion integrals for symmetric/asymmetric root cases.
//! * [`phitime`]  - coordinate time and azimuth (GEOPHITIME and helpers).
//! * [`radial`]   - GEOR: solve for the final radius given the final mu.
//! * [`geomu`]    - GEOMU: solve for the final mu given the final radius,
//!   including root classification (NCASE) and the M(mu) integrals.
//! * [`camera`]   - INITIALIZE_CAMERA_GEOKERR: camera pixel -> impact
//!   parameters and constants of motion.
//! * this file    - GEOKERR: the main per-geodesic driver.
//!
//! All routines are scalar; the Fortran versions were vectorized over a
//! single geodesic but every operation is element-wise, so scalar Rust is
//! numerically identical.

pub mod camera;
pub mod driver;
pub mod elliptic;
pub mod geomu;
pub mod imu;
pub mod phitime;
pub mod radial;
pub mod special;

pub use camera::{initialize_camera_geokerr, CameraArgs, CameraPixel};
pub use driver::{geokerr, GeokerrResult};
pub use geomu::{findmroots, indep_muf, GeomuResult};

/// Result of one geokerr call: arrays of the same length as 1 + npts + next - 2*kext
/// (matching upstream's `UFI, MUFI, DTI, DPHI, TPMI, TPRI, LAMBDAI`).
#[derive(Clone, Debug, Default)]
pub struct GeodesicPoints {
    /// inverse radius u = 1/r
    pub ufi: Vec<f64>,
    /// mu = cos(theta)
    pub mufi: Vec<f64>,
    /// delta t
    pub dti: Vec<f64>,
    /// delta phi
    pub dphi: Vec<f64>,
    /// number of mu turning points
    pub tpmi: Vec<i32>,
    /// number of u turning points
    pub tpri: Vec<i32>,
    /// affine parameter
    pub lambdai: Vec<f64>,
}
