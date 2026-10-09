#![allow(clippy::chunks_exact_to_as_chunks)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::type_complexity)]
#![allow(dead_code)]
//! Standard NT73 thin disk model.
//!
//! Direct translation of `fluid_model_thindisk.f90` (upstream GRTRANS).
//! All internal arithmetic is single precision (`real`) exactly as upstream;
//! values are widened to f64 only where upstream stores them in
//! `real(kind=8)` four-vector components.

use grtrans_core::constants::{C2, G, MSUN, PI, SIGB};
use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_cov_f32, calc_polvec, calc_rms, krolikc, ledd};

/// Module state (upstream `mbh, mdot, rin, rout`).
#[derive(Clone, Copy, Debug)]
pub struct ThindiskState {
    pub mbh: f32,
    pub mdot: f32,
    pub rin: f32,
    pub rout: f32,
}

impl Default for ThindiskState {
    fn default() -> Self {
        ThindiskState {
            mbh: 10.0,
            mdot: 0.1,
            rin: -1.0,
            rout: 1e8,
        }
    }
}

/// Upstream `init_thindisk`.
pub fn init_thindisk(state: &mut ThindiskState, a: f32, mdot: f32, mbh: f32, rin: f32, rout: f32) {
    let _ = a;
    state.mdot = mdot;
    state.mbh = mbh;
    state.rin = rin;
    state.rout = rout;
}

/// Upstream `thindisk_vals`: returns `(T, omega)` per input point.
///
/// `state.rin` is updated to `max(rms, rin)` as upstream does.
pub fn thindisk_vals(
    state: &mut ThindiskState,
    r: &[f32],
    th: &[f32],
    a: f32,
) -> (Vec<f32>, Vec<f32>) {
    let n = r.len();
    assert_eq!(th.len(), n);
    let mbh = state.mbh;
    let mdot = state.mdot;
    let msun_f = MSUN as f32; // f32 rounding of the f64 constant, as upstream
                              // lbh = MBH*MSUN*G/C2: mixed-kind expression evaluated in f64 then
                              // stored in a default real
    let lbh: f32 = ((mbh as f64) * MSUN * G / C2) as f32;
    let rms = calc_rms(a as f64) as f32;
    state.rin = rms.max(state.rin);
    let rin = state.rin;
    let rout = state.rout;
    let _ = msun_f;

    // kc = krolikc(r, a) (single precision)
    let mut t = vec![0.0f32; n];
    let mut omega = vec![0.0f32; n];
    let mdotedd: f32 = (ledd(mbh) / C2) as f32;
    // T0 = (3/8/pi*G*MBH*MSUN*Mdot*Mdotedd/lbh/lbh/lbh/SIGB)**0.25
    // (mixed-kind expression in f64, result stored in a default real)
    let t0_inner = 3.0f64 / 8.0 / PI * G * (mbh as f64) * MSUN * (mdot as f64) * (mdotedd as f64)
        / (lbh as f64)
        / (lbh as f64)
        / (lbh as f64)
        / SIGB;
    let t0: f32 = t0_inner.powf(0.25) as f32;

    for i in 0..n {
        let ri = r[i];
        let thi = th[i];
        let b = 1.0f32 - 3.0 / ri + 2.0 * a / ri.powf(3.0 / 2.0);
        let kc = krolikc(ri as f64, a as f64) as f32;
        let d = ri * ri - 2.0 * ri + a * a;
        let lc =
            (rms * rms - 2.0 * a * rms.sqrt() + a * a) / (rms.powf(1.5) - 2.0 * rms.sqrt() + a);
        let hc = (2.0 * ri - a * lc) / d;
        let ar = (ri * ri + a * a).powf(2.0) - a * a * d * thi.sin().powi(2);
        let om = 2.0 * a * ri / ar;
        omega[i] = if ri > rms {
            (1.0 / (ri.powf(3.0 / 2.0) + a)).max(om)
        } else {
            ((lc + a * hc) / (ri * ri + 2.0 * ri * (1.0 + hc))).max(om)
        };
        t[i] = if ri > rin && ri < rout {
            // T=T0*(kc/b/r**3.)**(1./4.) evaluated in single precision
            t0 * (kc / b / ri.powf(3.0)).powf(1.0 / 4.0)
        } else {
            (t0 as f64 / 1e5) as f32
        };
    }
    (t, omega)
}

/// Upstream `get_thindisk_fluidvars`.
pub fn get_thindisk_fluidvars(
    state: &mut ThindiskState,
    x0: &[FourVector],
    k0: &[FourVector],
    a: f32,
    f: &mut crate::fluid::Fluid,
) {
    let n = x0.len();
    let r32: Vec<f32> = x0.iter().map(|v| v.data[1] as f32).collect();
    let th32: Vec<f32> = x0.iter().map(|v| v.data[2] as f32).collect();
    let (t, omega) = thindisk_vals(state, &r32, &th32, a);
    f.rho = t;
    for i in 0..n {
        f.u[i].data[1] = 0.0;
        f.u[i].data[2] = 0.0;
        for c in 0..4 {
            f.b[i].data[c] = 0.0;
        }
        let metric = blmetric_cov_f32(r32[i], th32[i], a);
        let m = metric.0;
        f.u[i].assign_metric(metric);
        f.b[i].assign_metric(metric);
        let om = omega[i] as f64;
        let u1 = (-1.0 / (m[0] + 2.0 * m[3] * om + m[9] * om * om)).sqrt();
        f.u[i].data[0] = u1;
        f.u[i].data[3] = om * u1;
        // magnetic field proxy: polarization vector with psi = pi/2
        f.b[i] = calc_polvec(
            x0[i].data[1],
            x0[i].data[2].cos(),
            &k0[i],
            a as f64,
            std::f64::consts::FRAC_PI_2,
        );
    }
}
