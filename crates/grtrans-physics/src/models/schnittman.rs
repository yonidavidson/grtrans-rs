//! Hot spot model of Schnittman & Bertschinger (2004),
//! `fluid_model_hotspot_schnittman.f90` plus
//! `get_schnittman_hotspot_fluidvars` / `convert_fluidvars_schnittman_hotspot`
//! from `fluid.f90`. Cartesian Gaussian spot with a local Keplerian velocity
//! and a toroidal field (the f32 metric path, as upstream).

use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_cov_f32, calc_rms};

use crate::fluid::Fluid;

/// Upstream module state (`rspot, r0spot, n0spot, tspot`).
#[derive(Clone, Copy, Debug)]
pub struct SchnittmanState {
    pub rspot: f64,
    pub r0spot: f64,
    pub n0spot: f64,
    pub tspot: f64,
}

impl Default for SchnittmanState {
    fn default() -> Self {
        SchnittmanState {
            rspot: 1.5,
            r0spot: 6.0,
            n0spot: 1e4,
            tspot: 0.0,
        }
    }
}

/// Upstream `init_schnittman_hotspot`.
pub fn init_schnittman_hotspot(state: &mut SchnittmanState, rs: f64, r0: f64, n0: f64) {
    state.rspot = rs;
    state.r0spot = r0;
    state.n0spot = n0;
    state.tspot = 0.0;
}

/// Upstream `schnittman_hotspot_vals`: Cartesian Gaussian spot density.
pub fn schnittman_hotspot_vals(state: &SchnittmanState, x0: &[FourVector], a: f64) -> Vec<f64> {
    let omega = 1.0 / (state.r0spot.powf(1.5) + a);
    let dcrit = 16.0 * state.rspot * state.rspot;
    x0.iter()
        .map(|x| {
            let (r, phi, t, th) = (x.data[1], x.data[3], x.data[0], x.data[2]);
            let (sx, sy, sz) = (
                r * th.sin() * phi.cos(),
                r * th.sin() * phi.sin(),
                r * th.cos(),
            );
            let phispot = omega * (t + state.tspot);
            let (xspot, yspot) = (state.r0spot * phispot.cos(), state.r0spot * phispot.sin());
            let d2 = (sx - xspot).powi(2) + (sy - yspot).powi(2) + sz * sz;
            if d2 < dcrit {
                // upstream returns a unit Gaussian (n0spot is not used)
                (-d2 / (2.0 * state.rspot * state.rspot)).exp()
            } else {
                0.0
            }
        })
        .collect()
}

/// Upstream `get_schnittman_hotspot_fluidvars`: local Keplerian velocity with
/// the ISCO-branch merge, and a toroidal field (f32 metric).
pub fn get_schnittman_hotspot_fluidvars(
    state: &SchnittmanState,
    x0: &[FourVector],
    a: f64,
    f: &mut Fluid,
) {
    let n = schnittman_hotspot_vals(state, x0, a);
    let rms = calc_rms(a as f32 as f64) as f32 as f64;
    for (i, x) in x0.iter().enumerate() {
        let (r, th) = (x.data[1], x.data[2]);
        f.rho[i] = n[i] as f32;
        let omega = 1.0 / (r.powf(1.5) + a);
        let metric = blmetric_cov_f32(r as f32, th as f32, a as f32);
        let m = &metric.0;
        let safe = m[0] + 2.0 * m[3] * omega + m[9] * omega * omega;
        let d = r * r - 2.0 * r + a * a;
        let lc =
            (rms * rms - 2.0 * a * rms.sqrt() + a * a) / (rms.powf(1.5) - 2.0 * rms.sqrt() + a);
        let hc = (2.0 * r - a * lc) / d;
        let ar = (r * r + a * a).powi(2) - a * a * d * th.sin().powi(2);
        let om = 2.0 * a * r / ar;
        let omt = if r > rms {
            (1.0 / (r.powf(1.5) + a)).max(om)
        } else {
            ((lc + a * hc) / (r * r + 2.0 * r * (1.0 + hc))).max(om)
        };
        let ut = (-1.0 / (m[0] + 2.0 * m[3] * omt + m[9] * omt * omt)).sqrt();
        // upstream merge: safe<0 keeps the local-omega velocity, safe>=0
        // switches to the omt velocity
        let (u0, u3) = if safe < 0.0 {
            let u0 = (-1.0 / safe).sqrt();
            (u0, omega * u0)
        } else {
            (ut, omt * ut)
        };
        let gfac = 1.0
            / ((m[9] * m[0] - m[3] * m[3]) * (m[9] * u3 * u3 + u0 * (2.0 * m[3] * u3 + m[0] * u0)))
                .sqrt();
        let b0 = gfac * (m[9] * u3 + m[3] * u0).abs();
        let b3 = -(m[9] * u3 + m[3] * u0).signum() * (u0 * m[0] + m[3] * u3) * gfac;
        f.u[i] = FourVector::flat([u0, 0.0, 0.0, u3]);
        f.u[i].assign_metric(metric);
        f.b[i] = FourVector::flat([b0, 0.0, 0.0, b3]);
        f.b[i].assign_metric(metric);
        f.bmag[i] = f.b[i].dot(&f.b[i]).abs().sqrt() as f32;
    }
}

/// Upstream `convert_fluidvars_schnittman_hotspot`: ncgs=rho, bcgs=1.0.
/// Upstream leaves ncgsnth/tcgs unassigned; we set ncgsnth=rho and tcgs=0 so
/// the non-thermal emissivities are well-defined (documented choice).
pub fn convert_schnittman(f: &Fluid) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ncgs: Vec<f64> = f.rho.iter().map(|v| *v as f64).collect();
    let ncgsnth = ncgs.clone();
    let bcgs = vec![1.0; ncgs.len()];
    let tcgs = vec![0.0; ncgs.len()];
    (ncgs, ncgsnth, bcgs, tcgs)
}
