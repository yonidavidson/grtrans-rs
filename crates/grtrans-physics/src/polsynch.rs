//! Polarized synchrotron and blackbody emissivity helpers.
//!
//! Translation target: `polsynchemis.f90`. This module currently provides
//! the Planck function and will grow the thermal/power-law synchrotron
//! coefficients.

use grtrans_core::constants::{C2, H, K};

/// Planck spectrum for an array of temperatures (cgs), upstream `bnu`
/// (`polsynchemis.f90` lines 1014-1032).
///
/// Uses the Rayleigh-Jeans limit for `h nu / K T < 1e-6` and replaces exact
/// zeros by machine epsilon, as upstream.
pub fn bnu(t: &[f64], nu: f64) -> Vec<f64> {
    let mut out: Vec<f64> = t
        .iter()
        .map(|&temp| {
            if H * nu / K / temp < 1e-6 {
                2.0 * nu * nu * K * temp / C2
            } else {
                2.0 * H * nu / C2 * nu * nu / ((H * nu / K / temp).exp() - 1.0)
            }
        })
        .collect();
    for v in out.iter_mut() {
        if *v == 0.0 {
            *v = f64::EPSILON;
        }
    }
    out
}
