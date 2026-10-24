//! Semi-analytic RIAF model (Yuan, Quataert & Narayan 2003; Broderick+ 2009).
//!
//! Direct translation of `fluid_model_sariaf.f90` and the
//! `get_sariaf_fluidvars` / `convert_fluidvars_sariaf` routines from
//! `fluid.f90`. Single precision internally as upstream (with the same
//! mixed-kind narrowing points).

use grtrans_core::constants::{C2, MP, PI};
use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_cov_f32, calc_rms, calc_u0, rms_vel};

/// Upstream module state (`riaf_*`).
#[derive(Clone, Copy, Debug)]
pub struct SariafState {
    pub pnth: f32,
    pub n0: f32,
    pub t0: f32,
    pub nnth0: f32,
    pub beta: f32,
    pub bl06: i32,
}

impl Default for SariafState {
    fn default() -> Self {
        SariafState {
            pnth: 2.9,
            n0: 3e7,
            t0: 1.7e11,
            nnth0: 8e4,
            beta: 10.0,
            bl06: 0,
        }
    }
}

/// Upstream `init_sariaf` (optional args there; the caller supplies them).
pub fn init_sariaf(
    state: &mut SariafState,
    n0: f32,
    t0: f32,
    nnth0: f32,
    pnth: f32,
    beta: f32,
    bl06: i32,
) {
    state.n0 = n0;
    state.t0 = t0;
    state.nnth0 = nnth0;
    state.pnth = pnth;
    state.beta = beta;
    state.bl06 = bl06;
}

/// Upstream `sariaf_vals`: returns `(neth, te, B, vr, vth, omega, nenth)`.
#[allow(clippy::type_complexity)]
pub fn sariaf_vals(
    state: &SariafState,
    a: f32,
    u: &[f32],
    mu: &[f32],
) -> (
    Vec<f32>,
    Vec<f32>,
    Vec<f32>,
    Vec<f32>,
    Vec<f32>,
    Vec<f32>,
    Vec<f32>,
) {
    let n = u.len();
    let mut neth = vec![0.0f32; n];
    let mut nenth = vec![0.0f32; n];
    let mut te = vec![0.0f32; n];
    let mut b = vec![0.0f32; n];
    let mut vr = vec![0.0f32; n];
    let mut vth = vec![0.0f32; n];
    let mut omega = vec![0.0f32; n];
    for i in 0..n {
        let r = 1.0f32 / u[i];
        let rs = r / 2.0;
        let z = r * mu[i];
        let a2 = (r * r - z * z).sqrt();
        if state.bl06 != 1 {
            neth[i] = state.n0 * rs.powf(-1.1) * (-0.5 * (z / a2).powi(2)).exp();
            nenth[i] = state.nnth0 * rs.powf(-state.pnth) * (-0.5 * (z / a2).powi(2)).exp();
            te[i] = state.t0 * rs.powf(-0.84);
        } else {
            neth[i] = state.n0 * a2.powf(-1.1) * (-0.5 * (z / a2).powi(2)).exp();
            nenth[i] = state.nnth0 * a2.powf(-state.pnth) * (-0.5 * (z / a2).powi(2)).exp();
            te[i] = state.t0 * r.powf(-0.84);
        }
        b[i] = ((8.0 * PI * (neth[i] as f64) * MP * C2 / (rs as f64) / 12.0 / (state.beta as f64))
            .sqrt()) as f32;
        vr[i] = 0.0 * r;
        vth[i] = 0.0 * r;
        omega[i] = 1.0 / (r.powf(3.0 / 2.0) + a);
    }
    (neth, te, b, vr, vth, omega, nenth)
}

/// Upstream `get_sariaf_fluidvars`.
pub fn get_sariaf_fluidvars(
    state: &SariafState,
    x0: &[FourVector],
    a: f32,
    f: &mut crate::fluid::Fluid,
) {
    let n = x0.len();
    let rr: Vec<f32> = x0.iter().map(|v| v.data[1] as f32).collect();
    let ctheta: Vec<f32> = x0.iter().map(|v| v.data[2].cos() as f32).collect();
    let u: Vec<f32> = rr.iter().map(|v| 1.0 / v).collect();
    let (nn, tt, bmag, vr, vth, omega, nnth) = sariaf_vals(state, a, &u, &ctheta);
    let rms = calc_rms(a as f64) as f32;
    for i in 0..n {
        // upstream: kerr_metric(real(rr), real(th), a) -> single precision
        let metric = blmetric_cov_f32(rr[i], x0[i].data[2] as f32, a);
        let u0 = calc_u0(&metric, vr[i] as f64, vth[i] as f64, omega[i] as f64);
        let fu = rms_vel(a as f64, x0[i].data[2], x0[i].data[1]);
        if rr[i] < rms {
            f.u[i].data[0] = fu.data[0];
            f.u[i].data[1] = fu.data[1];
            f.u[i].data[2] = fu.data[2];
            f.u[i].data[3] = fu.data[3];
        } else {
            f.u[i].data[0] = u0;
            f.u[i].data[1] = vr[i] as f64 * u0;
            f.u[i].data[2] = vth[i] as f64 * u0;
            f.u[i].data[3] = omega[i] as f64 * u0;
        }
        let m = &metric.0;
        let aleph = (-(m[3] * f.u[i].data[0] + m[9] * f.u[i].data[3])
            / (m[0] * f.u[i].data[0] + m[3] * f.u[i].data[3])) as f32;
        let bb =
            (m[0] * (aleph as f64) * (aleph as f64) + m[9] + 2.0 * m[3] * (aleph as f64)) as f32;
        let b4 = bmag[i] / bb.sqrt();
        f.b[i].data[0] = aleph as f64 * b4 as f64;
        f.b[i].data[1] = 0.0;
        f.b[i].data[2] = 0.0;
        f.b[i].data[3] = b4 as f64;
        f.u[i].assign_metric(metric);
        f.b[i].assign_metric(metric);
    }
    f.rho = nn;
    f.p = tt;
    f.bmag = bmag;
    f.rho2 = nnth;
}
