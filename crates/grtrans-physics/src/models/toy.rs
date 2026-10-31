//! Toy fluid model for the code-comparison tests (`fluid_model_toy.f90`)
//! plus `get_toy_fluidvars` / `convert_fluidvars_toy` from `fluid.f90`.
//!
//! Upstream defects (see docs/VALIDATION_PLAN.md §3.6): `aleph` is computed
//! from the *previous* call's `f%u` (stale/uninitialized) and `bb` is never
//! assigned, so the magnetic field is garbage. We make two documented
//! deterministic choices to obtain a well-defined model: compute the
//! four-velocity first and evaluate `aleph` from it, and use the standard
//! toroidal normalization `bb = g_tt*aleph^2 + g_phiphi + 2*g_tphi*aleph`
//! (as in SARIAF/POWERLAW). The density profile and the four-velocity
//! construction are well-defined upstream and are fixture-validated.

use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_con, blmetric_cov};

use crate::fluid::Fluid;

/// Upstream module state (`toy_l0, toy_n0, toy_h, toy_which`).
#[derive(Clone, Copy, Debug)]
pub struct ToyState {
    pub l0: f64,
    pub n0: f64,
    pub h: f64,
    pub which: i32,
}

impl Default for ToyState {
    fn default() -> Self {
        ToyState {
            l0: 10.0,
            n0: 3e7,
            h: 1.7,
            which: 0,
        }
    }
}

/// Upstream `init_toy` (n0, h, l0, which).
pub fn init_toy(state: &mut ToyState, n0: f64, h: f64, l0: f64, which: i32) {
    state.n0 = n0;
    state.h = h;
    state.l0 = l0;
    state.which = which;
}

/// Upstream `toy_vals` (code-comparison equations 1 and 2): density `n` and
/// specific angular momentum `l`.
pub fn toy_vals(state: &ToyState, a: f64, mu: &[f64], r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let q = 0.5;
    let _ = a;
    mu.iter()
        .zip(r.iter())
        .map(|(&m, &ri)| {
            let rcyl = ri * (1.0 - m * m).sqrt();
            let l = state.l0 / (1.0 + rcyl) * rcyl.powf(1.0 + q);
            let z = state.h * m;
            let dist = (ri / 10.0).powi(2) + z * z;
            let n = if dist < 20.0 {
                state.n0 * (-dist / 2.0).exp()
            } else {
                0.0
            };
            (n, l)
        })
        .unzip()
}

/// Upstream `get_toy_fluidvars` with the documented aleph/bb choices (see the
/// module doc comment).
pub fn get_toy_fluidvars(state: &ToyState, x0: &[FourVector], a: f64, f: &mut Fluid) {
    let rr: Vec<f64> = x0.iter().map(|v| v.data[1]).collect();
    let ctheta: Vec<f64> = x0.iter().map(|v| v.data[2].cos()).collect();
    let (n, l) = toy_vals(state, a, &ctheta, &rr);
    for (i, x) in x0.iter().enumerate() {
        let (r, th) = (x.data[1], x.data[2]);
        let metric = blmetric_cov(r, th, a);
        let con = blmetric_con(r, th, a);
        let (m, cm) = (&metric.0, &con.0);
        // four-velocity from the contravariant construction (well-defined
        // upstream): ubar/ulower, then f%u. Compute it BEFORE aleph.
        let li = l[i];
        let ubar = (-1.0 / (cm[0] + li * li * cm[9] - 2.0 * li * cm[3])).sqrt();
        let (ul_t, ul_phi) = (-ubar, li * ubar);
        let u0 = cm[3] * ul_phi + cm[0] * ul_t;
        let u3 = cm[3] * ul_t + cm[9] * ul_phi;
        // documented choice: aleph from the computed four-velocity
        let aleph = -(m[3] * u0 + m[9] * u3) / (m[0] * u0 + m[3] * u3);
        // documented choice: standard toroidal normalization for bb
        let bb = m[0] * aleph * aleph + m[9] + 2.0 * m[3] * aleph;
        let bmag = 1.0;
        let b4 = bmag / bb.sqrt();
        f.u[i] = FourVector::flat([u0, 0.0, 0.0, u3]);
        f.u[i].assign_metric(metric);
        f.b[i] = FourVector::flat([aleph * b4, 0.0, 0.0, b4]);
        f.b[i].assign_metric(metric);
        f.rho[i] = n[i] as f32;
        f.bmag[i] = bmag as f32;
    }
}

/// Upstream `convert_fluidvars_toy`: ncgs=rho, bcgs=bmag, tcgs=0. We set
/// ncgsnth=rho as well so non-thermal emissivities are well-defined
/// (documented choice; upstream leaves it unassigned).
pub fn convert_toy(f: &Fluid) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ncgs: Vec<f64> = f.rho.iter().map(|v| *v as f64).collect();
    let ncgsnth = ncgs.clone();
    let bcgs: Vec<f64> = f.bmag.iter().map(|v| *v as f64).collect();
    let tcgs = vec![0.0; ncgs.len()];
    (ncgs, ncgsnth, bcgs, tcgs)
}
