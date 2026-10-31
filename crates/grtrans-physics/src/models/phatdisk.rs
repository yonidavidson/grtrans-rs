//! Novikov–Thorne thin disk with a multicolor (log-normal broadened) blackbody
//! frequency table (`fluid_model_phatdisk.f90`), plus the
//! `get_phatdisk_fluidvars` / `convert_fluidvars_phatdisk` routines from
//! `fluid.f90`. The per-point flux table feeds the INTERP emissivity.
//!
//! Single precision internally as upstream, with the same mixed-kind
//! narrowing points. `r_tab` deliberately starts at index `i=1` (not `i-1`),
//! an upstream off-by-one that shifts the radius grid.

use grtrans_core::constants::{C2, H, K as KB, PI};
use grtrans_core::four_vector::FourVector;
use grtrans_core::interpolate::get_weight_f32;
use grtrans_core::kerr::{blmetric_cov_f32, calc_polvec};
use grtrans_core::math::tsum;

use crate::fluid::Fluid;
use crate::models::thindisk::{thindisk_vals, ThindiskState};

/// Upstream module state (namelist parameters + computed tables).
#[derive(Clone, Debug)]
pub struct PhatdiskState {
    pub nw: i32,
    pub wmin: f32,
    pub wmax: f32,
    pub nfreq_tab: i32,
    pub fmin: f32,
    pub fmax: f32,
    pub rmax: f32,
    pub nr: i32,
    pub sigt: f32,
    pub fcol: f32,
    /// thindisk parameters (MBH, Mdot, rin, rout). Upstream reads these from a
    /// thindisk namelist file; we take them from `FluidArgs` (documented
    /// choice, with the same values).
    pub thin: ThindiskState,
    // tables (filled by `init_phatdisk`)
    pub w: Vec<f32>,
    pub freq_tab: Vec<f32>,
    pub r_tab: Vec<f32>,
    pub om_tab: Vec<f32>,
    /// `fnu_tab(i,k)` flattened as `i*nfreq_tab + k` (radius-major, freq-fast)
    pub fnu_tab: Vec<f32>,
}

impl PhatdiskState {
    pub fn new(thin: ThindiskState) -> Self {
        PhatdiskState {
            nw: 500,
            wmin: 1e-4,
            wmax: 1e4,
            nfreq_tab: 100,
            fmin: 1e10,
            fmax: 1e15,
            rmax: 1e4,
            nr: 500,
            sigt: 0.4,
            fcol: 1.7,
            thin,
            w: Vec::new(),
            freq_tab: Vec::new(),
            r_tab: Vec::new(),
            om_tab: Vec::new(),
            fnu_tab: Vec::new(),
        }
    }
}

/// Upstream `init_phatdisk`: build the frequency/radius tables and the
/// multicolor blackbody flux table.
pub fn init_phatdisk(state: &mut PhatdiskState, a: f32) {
    let nw = state.nw as usize;
    let nft = state.nfreq_tab as usize;
    let nr = state.nr as usize;
    // w, freq_tab log-spaced (i-1 convention)
    state.w = (0..nw)
        .map(|i| {
            if nw == 1 {
                state.wmin
            } else {
                state.wmin * ((i as f32) / ((nw - 1) as f32) * (state.wmax / state.wmin).ln()).exp()
            }
        })
        .collect();
    state.freq_tab = (0..nft)
        .map(|i| {
            if nft == 1 {
                state.fmin
            } else {
                state.fmin
                    * ((i as f32) / ((nft - 1) as f32) * (state.fmax / state.fmin).ln()).exp()
            }
        })
        .collect();
    // r_tab: upstream uses (i, i=1, nr), i.e. starting at 1 (off-by-one)
    let rh = 1.0 + (1.0 - a * a).sqrt();
    state.r_tab = (0..nr)
        .map(|i| {
            let ii = (i + 1) as f32; // 1-based i
            rh * (ii / ((nr - 1) as f32) * (state.rmax / rh).ln()).exp()
        })
        .collect();
    // thin-disk temperature and omega on r_tab (theta = 0)
    let zeros = vec![0.0f32; nr];
    let (t, om) = thindisk_vals(&mut state.thin, &state.r_tab, &zeros, a);
    state.om_tab = om;
    // log-normal flux weights
    let x: Vec<f32> = state.w.iter().map(|v| v.ln()).collect();
    let l10sigt = std::f32::consts::LN_10 * state.sigt;
    let fw: Vec<f32> = x
        .iter()
        .map(|xi| {
            let arg = (xi + l10sigt * l10sigt) / l10sigt;
            (-arg * arg).exp() / l10sigt / PI.sqrt() as f32
        })
        .collect();
    state.fnu_tab = vec![0.0f32; nr * nft];
    for k in 0..nft {
        for i in 0..nr {
            // z = h*freq_tab(k)/(kb*T*fcol): f64 expression narrowed to f32
            let z = ((H * (state.freq_tab[k] as f64)) / (KB * (t[i] as f64) * (state.fcol as f64)))
                as f32;
            let denom: Vec<f32> = (0..nw)
                .map(|j| {
                    let ratio = z / state.w[j];
                    if ratio > 1e-4 {
                        ratio.exp() - 1.0
                    } else {
                        ratio
                    }
                })
                .collect();
            let igrand: Vec<f32> = (0..nw)
                .map(|j| if fw[j] > 0.0 { fw[j] / denom[j] } else { 0.0 })
                .collect();
            let xf: Vec<f64> = x.iter().map(|v| *v as f64).collect();
            let yf: Vec<f64> = igrand.iter().map(|v| *v as f64).collect();
            let integral = tsum(&xf, &yf)[nw - 1];
            // fnu_tab(i,k): f64 expression narrowed to f32
            let val = (state.fcol as f64).powi(-4)
                * 2.0
                * PI
                * (z as f64).powi(3)
                * (KB * (state.fcol as f64) * (t[i] as f64)).powi(3)
                / H
                / H
                / C2
                * integral;
            state.fnu_tab[i * nft + k] = val as f32;
        }
    }
}

/// Upstream `phatdisk_vals`: interpolate the flux table and omega onto the
/// given radii. Returns `(fnu, omega)` with `fnu` flattened point-major,
/// freq-fastest (the layout the INTERP emissivity expects after upstream's
/// `reshape(transpose(fnu), ...)`).
pub fn phatdisk_vals(state: &PhatdiskState, r: &[f32], _a: f32) -> (Vec<f32>, Vec<f32>) {
    let nft = state.nfreq_tab as usize;
    let np = r.len();
    let mut fnu = vec![0.0f32; np * nft];
    let mut omega = vec![0.0f32; np];
    for (p, &ri) in r.iter().enumerate() {
        let (w, j) = get_weight_f32(&state.r_tab, ri, 0);
        let lo = (j - 1) as usize; // 0-based lower bracket
        let hi = j as usize; // 0-based upper bracket
        for k in 0..nft {
            fnu[p * nft + k] =
                (1.0 - w) * state.fnu_tab[lo * nft + k] + w * state.fnu_tab[hi * nft + k];
        }
        omega[p] = (1.0 - w) * state.om_tab[lo] + w * state.om_tab[hi];
    }
    (fnu, omega)
}

/// Upstream `get_phatdisk_fluidvars`: Keplerian velocity, the per-point flux
/// table and the polarization vector.
pub fn get_phatdisk_fluidvars(
    state: &PhatdiskState,
    x0: &[FourVector],
    k0: &[FourVector],
    a: f32,
    f: &mut Fluid,
) {
    let rr: Vec<f32> = x0.iter().map(|v| v.data[1] as f32).collect();
    let (fnu, omega) = phatdisk_vals(state, &rr, a);
    f.fnu = fnu;
    f.freq_tab = state.freq_tab.clone();
    f.nfreq = state.nfreq_tab as usize;
    for (i, x) in x0.iter().enumerate() {
        let th = x.data[2] as f32;
        let metric = blmetric_cov_f32(rr[i], th, a);
        let m = &metric.0;
        let om = omega[i] as f64;
        let u0 = (-1.0 / (m[0] + 2.0 * m[3] * om + m[9] * om * om)).sqrt();
        f.u[i].data[0] = u0;
        f.u[i].data[1] = 0.0;
        f.u[i].data[2] = 0.0;
        f.u[i].data[3] = om * u0;
        f.u[i].assign_metric(metric);
        f.b[i] = calc_polvec(x.data[1], x.data[2].cos(), &k0[i], a as f64, 0.0);
        f.b[i].assign_metric(metric);
    }
}
