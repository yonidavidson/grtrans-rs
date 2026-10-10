//! Bondi/Michel spherical accretion (GR version).
//!
//! Direct translation of `fluid_model_sphacc.f90` (upstream GRTRANS). The
//! stored flow solution (temperature and velocity profiles) is embedded
//! from the upstream source; all arithmetic is single precision as upstream.

use crate::sphacc_tables::{RVALS, TVALS, UVALS, VVALS};
use grtrans_core::constants::{C2, MP, PI};
use grtrans_core::four_vector::FourVector;
use grtrans_core::interpolate::get_weight_f32;

/// Upstream `temper`: temperature profile interpolated at u = 1/r.
fn temper(u: &[f32]) -> Vec<f32> {
    let uv: Vec<f32> = RVALS.iter().map(|v| 1.0 / (*v as f32)).collect();
    let tv: Vec<f32> = TVALS.iter().map(|v| *v as f32).collect();
    let mut lindx: isize = 1;
    let mut out = vec![0.0f32; u.len()];
    for (i, &ui) in u.iter().enumerate() {
        let (weight, j) = get_weight_f32(&uv, ui, lindx);
        lindx = j;
        let j0 = (j - 1) as usize;
        out[i] = (1.0 - weight) * tv[j0] + weight * tv[j0 + 1];
    }
    out
}

/// Upstream `vel`: velocity profile interpolated on the stored u grid.
fn vel(u: &[f32]) -> Vec<f32> {
    let uv: Vec<f32> = UVALS.iter().map(|v| *v as f32).collect();
    let vv: Vec<f32> = VVALS.iter().map(|v| *v as f32).collect();
    let mut lindx: isize = 1;
    let mut out = vec![0.0f32; u.len()];
    for (i, &ui) in u.iter().enumerate() {
        let (weight, j) = get_weight_f32(&uv, ui, lindx);
        lindx = j;
        let j0 = (j - 1) as usize;
        out[i] = (1.0 - weight) * vv[j0] + weight * vv[j0 + 1];
    }
    out
}

/// Upstream `sphacc_vals`: returns `(n, B, T, ur)`.
pub fn sphacc_vals(u: &[f32]) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    let gamma = 5.0f32 / 3.0;
    let alpha = 1.0f32 / 4.0;
    let ninf = 1.0f32;
    let tinf = 0.917e-9f32;
    let us = 0.94f32;
    let n: Vec<f32> = u
        .iter()
        .map(|&ui| ninf * alpha / 4.0 / us * (2.0 * gamma * tinf).powf(-1.5) * (2.0 * ui).powf(1.5))
        .collect();
    let t = temper(u);
    let b: Vec<f32> = (0..u.len())
        .map(|i| {
            // B=sqrt(8.*pi*n*mp/2.*c2*u): mixed-kind expression evaluated in
            // f64 and stored in a default real
            ((8.0 * PI * (n[i] as f64) * MP / 2.0 * C2 * (u[i] as f64)).sqrt()) as f32
        })
        .collect();
    let ur = vel(u);
    (n, b, t, ur)
}

/// Upstream `get_sphacc_fluidvars`.
pub fn get_sphacc_fluidvars(x0: &[FourVector], f: &mut crate::fluid::Fluid) {
    let n = x0.len();
    let u32: Vec<f32> = x0.iter().map(|v| 1.0 / v.data[1] as f32).collect();
    let (nn, b, t, ur) = sphacc_vals(&u32);
    f.rho = nn;
    f.p = t;
    for i in 0..n {
        let ui = u32[i];
        let g00 = -(1.0f32 - 2.0 * ui);
        let grr: f32 = (-1.0f64 / g00 as f64) as f32;
        f.u[i].data[1] = -(ur[i] as f64);
        f.u[i].data[2] = 0.0;
        f.u[i].data[3] = 0.0;
        let ur64 = f.u[i].data[1];
        let grr64 = grr as f64;
        let g0064 = g00 as f64;
        f.u[i].data[0] = ((-grr64 * ur64 * ur64 - 1.0) / g0064).sqrt();
        f.b[i].data[2] = 0.0;
        f.b[i].data[3] = 0.0;
        let b64 = b[i] as f64;
        let u0 = f.u[i].data[0];
        f.b[i].data[0] = (ur64 * ur64 * grr64 * b64 * b64
            / (ur64 * ur64 * g0064 * grr64 + u0 * u0 * g0064 * g0064))
            .sqrt();
        f.b[i].data[1] =
            -(b64 * b64 / grr64 - f.b[i].data[0] * f.b[i].data[0] * g0064 / grr64).sqrt();
    }
    f.bmag = b;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn velocity_profile_is_smooth() {
        // u values along a nearly radial SPHACC ray (r from ~200 to ~2)
        let u: Vec<f32> = (0..100)
            .map(|i| {
                let r = 200.50125f32 - (i as f32) * (198.49f32 / 99.0);
                1.0 / r
            })
            .collect();
        let ur = vel(&u);
        for (i, v) in ur.iter().enumerate() {
            assert!(
                v.is_finite() && v.abs() < 1.0,
                "vel[{i}] = {v} (u = {})",
                u[i]
            );
        }
        // the stored velocity profile increases monotonically with u
        for i in 1..ur.len() {
            assert!(
                ur[i] >= ur[i - 1] - 1e-6,
                "not monotone at {i}: {} -> {}",
                ur[i - 1],
                ur[i]
            );
        }
    }
}
