//! Gaussian hot spot on a Keplerian orbit (Broderick & Loeb 2005, 2006),
//! `fluid_model_hotspot.f90` plus `get_hotspot_fluidvars` /
//! `convert_fluidvars_hotspot` from `fluid.f90`.
//!
//! All quantities are computed in double precision (the upstream uses the f64
//! `kerr_metric` here, unlike the f32 path of SARIAF/PHATDISK). `tspot` is
//! initialized to zero (single-snapshot use, as in the regression setups).

use grtrans_core::constants::PI;
use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_cov, calc_polvec, calc_rms};

use crate::fluid::Fluid;

/// Upstream module state (`rspot, r0spot, n0spot, tspot, bl06`).
#[derive(Clone, Copy, Debug)]
pub struct HotspotState {
    pub rspot: f64,
    pub r0spot: f64,
    pub n0spot: f64,
    pub tspot: f64,
    pub bl06: i32,
}

impl Default for HotspotState {
    fn default() -> Self {
        HotspotState {
            rspot: 1.5,
            r0spot: 6.0,
            n0spot: 1e4,
            tspot: 0.0,
            bl06: 0,
        }
    }
}

/// Upstream `init_hotspot` (optional-arg form used by the driver).
pub fn init_hotspot(state: &mut HotspotState, rs: f64, r0: f64, n0: f64, bl: i32) {
    state.rspot = rs;
    state.r0spot = r0;
    state.n0spot = n0;
    state.bl06 = bl;
    state.tspot = 0.0;
}

/// Upstream `hotspot_vals`: returns `(n, b, u, xout)` where `n` is the spot
/// density, `b`/`u` the magnetic field and four-velocity, and `xout` the
/// spot-frame coordinates.
pub fn hotspot_vals(
    state: &HotspotState,
    x0: &[FourVector],
    k0: &[FourVector],
    a: f64,
) -> (Vec<f64>, Vec<FourVector>, Vec<FourVector>, Vec<FourVector>) {
    let np = x0.len();
    let rms = calc_rms(a as f32 as f64) as f32 as f64;
    let dcrit = 8.0;
    let bmin = 1e-6;
    let omega = 1.0 / (state.r0spot.powf(1.5) + a);
    // spot position and metric
    let xspot = FourVector::flat([0.0, state.r0spot, 0.0f64.acos(), 0.0]);
    let metrics = blmetric_cov(state.r0spot, xspot.data[2], a);
    let uspot_ut =
        (-1.0 / (metrics.0[0] + 2.0 * metrics.0[3] * omega + metrics.0[9] * omega * omega)).sqrt();
    let mut n = vec![0.0f64; np];
    let mut b = vec![FourVector::flat([0.0; 4]); np];
    let mut u = vec![FourVector::flat([0.0; 4]); np];
    let mut xout = Vec::with_capacity(np);
    for (i, xi) in x0.iter().enumerate() {
        let mut x = *xi;
        let r = x.data[1];
        let th = x.data[2];
        // rotate to the spot frame (x%data(4) already shifted by the caller)
        x.data[3] -= (state.tspot + x.data[0]) * omega;
        x.data[3] = x.data[3].sin().atan2(x.data[3].cos());
        let metric = blmetric_cov(r, th, a);
        let m = &metric.0;
        // separation four-vector dx = xspot - x
        let mut dx = FourVector::flat([
            0.0,
            xspot.data[1] - x.data[1],
            xspot.data[2] - x.data[2],
            xspot.data[3] - x.data[3],
        ]);
        dx.data[0] = 0.0;
        dx.assign_metric(metrics);
        // uspot at the spot (contravariant Keplerian velocity)
        let mut uspot = FourVector::flat([uspot_ut, 0.0, 0.0, omega * uspot_ut]);
        uspot.assign_metric(metrics);
        let dnorm = dx.dot(&dx) + dx.dot(&uspot).powi(2);
        // four-velocity (Keplerian with the ISCO branch and the safe guard)
        let d = r * r - 2.0 * r + a * a;
        let lc =
            (rms * rms - 2.0 * a * rms.sqrt() + a * a) / (rms.powf(1.5) - 2.0 * rms.sqrt() + a);
        let hc = (2.0 * r - a * lc) / d;
        let ar = (r * r + a * a).powi(2) - a * a * d * th.sin().powi(2);
        let om = 2.0 * a * r / ar;
        let safe = m[0] + 2.0 * m[3] * omega + m[9] * omega * omega;
        let omt = if r > rms {
            (1.0 / (r.powf(1.5) + a)).max(om)
        } else {
            ((lc + a * hc) / (r * r + 2.0 * r * (1.0 + hc))).max(om)
        };
        let ut = (-1.0 / (m[0] + 2.0 * m[3] * omt + m[9] * omt * omt)).sqrt();
        // upstream merge: safe<0 keeps the spot's constant-omega velocity,
        // safe>=0 switches to the local omt velocity
        let (u0, u3) = if safe < 0.0 {
            let u0 = (-1.0 / safe).sqrt();
            (u0, omega * u0)
        } else {
            (ut, omt * ut)
        };
        let _ = bmin;
        // Gaussian density (cut off beyond dcrit)
        let arg = dnorm / (2.0 * state.rspot * state.rspot);
        n[i] = if arg < dcrit {
            state.n0spot * (-arg).exp()
        } else {
            0.0
        };
        // equipartition field strength
        let bmag = if arg < dcrit {
            (0.1 * 8.0 * PI * n[i] * 100.0 * 1.67e-24 / 2.0 * 9e20 / r).sqrt()
        } else {
            1.0
        };
        u[i] = FourVector::flat([u0, 0.0, 0.0, u3]);
        u[i].assign_metric(metric);
        // magnetic field model (bl06: 1 toroidal, 0 poloidal, 2 vertical, 3 polvec)
        let bl = state.bl06;
        let mut bv = if bl.abs() == 1 {
            let gfac = 1.0
                / ((m[9] * m[0] - m[3] * m[3])
                    * (m[9] * u3 * u3 + u0 * (2.0 * m[3] * u3 + m[0] * u0)))
                    .sqrt();
            let b1 = bmag * gfac * (m[9] * u3 + m[3] * u0).abs();
            let b4 = -bmag * (m[9] * u3 + m[3] * u0).signum() * (u0 * m[0] + m[3] * u3) * gfac;
            FourVector::flat([b1, 0.0, 0.0, b4])
        } else if bl == 0 {
            FourVector::flat([0.0, 0.0, bmag / m[7].sqrt(), 0.0])
        } else if bl.abs() == 2 {
            FourVector::flat([
                0.0,
                -bmag / m[4].sqrt() * th.cos(),
                bmag / m[7].sqrt() * th.sin(),
                0.0,
            ])
        } else if bl.abs() == 3 {
            calc_polvec(r, th.cos(), &k0[i], a, 1.0f64.asin())
        } else {
            FourVector::flat([0.0; 4])
        };
        bv.assign_metric(metric);
        b[i] = bv;
        xout.push(x);
    }
    (n, b, u, xout)
}

/// Upstream `get_hotspot_fluidvars`: transforms phi/t, then `hotspot_vals`.
pub fn get_hotspot_fluidvars(
    state: &HotspotState,
    x0: &[FourVector],
    k0: &[FourVector],
    a: f64,
    f: &mut Fluid,
) {
    // upstream: x%data(4) = -1*acos(0.) - x%data(4); x%data(1) = -1.*x%data(1)
    let xs: Vec<FourVector> = x0
        .iter()
        .map(|xi| {
            let mut x = *xi;
            x.data[3] = -0.0f64.acos() - x.data[3];
            x.data[0] = -x.data[0];
            x
        })
        .collect();
    let (n, b, u, _xout) = hotspot_vals(state, &xs, k0, a);
    f.rho = n.iter().map(|v| *v as f32).collect();
    for i in 0..f.u.len() {
        f.u[i] = u[i];
        f.b[i] = b[i];
        f.bmag[i] = b[i].dot(&b[i]).abs().sqrt() as f32;
    }
}

/// Upstream `convert_fluidvars_hotspot`: ncgs=rho, bcgs=bmag, ncgsnth=rho.
pub fn convert_hotspot(f: &Fluid) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ncgs: Vec<f64> = f.rho.iter().map(|v| *v as f64).collect();
    let bcgs: Vec<f64> = f.bmag.iter().map(|v| *v as f64).collect();
    let ncgsnth = ncgs.clone();
    // upstream leaves tcgs unassigned (non-thermal emissivities only)
    let tcgs = vec![0.0; ncgs.len()];
    (ncgs, ncgsnth, bcgs, tcgs)
}
