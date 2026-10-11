//! Power-law jet model with a toroidal magnetic field.
//!
//! Direct translation of `fluid_model_powerlaw.f90` (upstream GRTRANS) and
//! the `get_powerlaw_fluidvars` / `convert_fluidvars_powerlaw` /
//! `monika_e` routines from `fluid.f90`. The model is single precision
//! internally except where upstream mixes kinds; the port reproduces the
//! same narrowing points.

use grtrans_core::constants::{C2, MP, PI};
use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_cov, calc_u0};

/// Upstream module state (`pl_*`).
#[derive(Clone, Copy, Debug)]
pub struct PowerlawState {
    pub pnth: f32,
    pub n0: f32,
    pub t0: f32,
    pub nnth0: f32,
    pub beta: f32,
    pub pn: f32,
    pub pt: f32,
    pub rin: f32,
    pub rout: f32,
    pub thin: f32,
    pub thout: f32,
    pub phiin: f32,
    pub phiout: f32,
}

impl Default for PowerlawState {
    fn default() -> Self {
        PowerlawState {
            pnth: 0.0,
            n0: 3e7,
            t0: 6e10,
            nnth0: 8e4,
            beta: 10.0,
            pn: 0.0,
            pt: 0.0,
            rin: 0.0,
            rout: 1e8,
            thin: -10.0,
            thout: 10.0,
            phiin: 0.0,
            phiout: 1e4,
        }
    }
}

/// Upstream `init_powerlaw` (all arguments optional there).
#[allow(clippy::too_many_arguments)]
pub fn init_powerlaw(
    state: &mut PowerlawState,
    n0: f32,
    t0: f32,
    nnth0: f32,
    pnth: f32,
    beta: f32,
    pn: f32,
    pt: f32,
    rin: f32,
    rout: f32,
    thin: f32,
    thout: f32,
    phiin: f32,
    phiout: f32,
) {
    state.n0 = n0;
    state.t0 = t0;
    state.nnth0 = nnth0;
    state.pnth = pnth;
    state.beta = beta;
    state.pn = pn;
    state.pt = pt;
    state.rin = rin;
    state.rout = rout;
    state.thin = thin;
    state.thout = thout;
    state.phiin = phiin;
    state.phiout = phiout;
}

/// Upstream `powerlaw_vals`: returns `(neth, te, B, vr, vth, omega, nenth)`.
#[allow(clippy::type_complexity)]
pub fn powerlaw_vals(
    state: &PowerlawState,
    u: &[f32],
    mu: &[f32],
) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    let n = u.len();
    let mut neth = vec![0.0f32; n];
    let mut nenth = vec![0.0f32; n];
    let mut te = vec![0.0f32; n];
    let mut vr = vec![0.0f32; n];
    let mut vth = vec![0.0f32; n];
    let mut omega = vec![0.0f32; n];
    let mut b = vec![0.0f32; n];
    for i in 0..n {
        let r = 1.0f32 / u[i];
        let rs = r / 2.0;
        neth[i] = state.n0 * rs.powf(-state.pn);
        nenth[i] = state.nnth0 * rs.powf(-state.pnth);
        te[i] = state.t0 * rs.powf(-state.pt);
        vr[i] = 0.0 * r;
        vth[i] = 0.0 * r;
        omega[i] = state.phiin / r;
        if r > state.rout || r < state.rin || mu[i] < state.thin || mu[i] > state.thout {
            neth[i] = 0.0;
            nenth[i] = 0.0;
        }
        // mixed-kind expression evaluated in f64 and stored in a default real
        b[i] = ((8.0 * PI * (neth[i] as f64) * MP * C2 / 10.0 / 12.0 / (state.beta as f64)).sqrt())
            as f32;
    }
    (neth, te, b, vr, vth, omega, nenth)
}

/// Upstream `monika_e`: Moscibrodzka+2016 e- model (`rlow`/`rhigh` are
/// T_p/T_e bounds).
pub fn monika_e(rho: &[f32], p: &[f32], b: &[f32], rlow: f64, rhigh: f64) -> Vec<f64> {
    let n = rho.len();
    let mut trat = vec![0.0f64; n];
    for i in 0..n {
        // beta = p/(b*b)/0.5d0: p/(b*b) in f32, then f64 division
        let beta = (p[i] / (b[i] * b[i])) as f64 / 0.5;
        let b2 = beta * beta;
        trat[i] = if b[i] > 0.0 {
            rhigh * b2 / (1.0 + b2) + rlow / (1.0 + b2)
        } else {
            rhigh
        };
    }
    let _ = rho;
    trat
}

/// Upstream `get_powerlaw_fluidvars`.
pub fn get_powerlaw_fluidvars(
    state: &PowerlawState,
    x0: &[FourVector],
    a: f32,
    f: &mut crate::fluid::Fluid,
) {
    let n = x0.len();
    // f32 arrays as upstream
    let rr: Vec<f32> = x0.iter().map(|v| v.data[1] as f32).collect();
    let ctheta: Vec<f32> = x0.iter().map(|v| v.data[2].cos() as f32).collect();
    let u: Vec<f32> = rr.iter().map(|v| 1.0 / v).collect();
    let (nn, tt, bmag, vr, vth, omega, nnth) = powerlaw_vals(state, &u, &ctheta);
    for i in 0..n {
        let metric = blmetric_cov(x0[i].data[1], x0[i].data[2], a as f64);
        let u0 = calc_u0(&metric, vr[i] as f64, vth[i] as f64, omega[i] as f64);
        f.u[i].data[0] = u0;
        f.u[i].data[1] = vr[i] as f64 * u0;
        f.u[i].data[2] = vth[i] as f64 * u0;
        f.u[i].data[3] = omega[i] as f64 * u0;
        // aleph and bb are default-real intermediates upstream (f64
        // expressions narrowed to f32)
        let m = &metric.0;
        let aleph = (-(m[3] * u0 + m[9] * f.u[i].data[3])
            / (m[0] * u0 + m[3] * f.u[i].data[3])) as f32;
        let bb = (m[0] * (aleph as f64) * (aleph as f64) + m[9] + 2.0 * m[3] * (aleph as f64)) as f32;
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
