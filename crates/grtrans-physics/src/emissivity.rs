//! Emissivity dispatch and blackbody models.
//!
//! Direct translation of `emis.f90` (upstream GRTRANS): emissivity types,
//! the `emis` object, the K-buffer (11-output) convention, blackbody
//! models (BB/FBB/BBPOL), and the frame rotations/scalings applied by the
//! driver.

use crate::polsynch::{bnu, polsynchpl, polsynchth, synchemis, synchpl};
use grtrans_core::chandra::interp_chandra_tab24_f32;

/// Upstream emissivity type constants (`emis.f90` lines 10-17).
pub mod etype {
    pub const ELAMBDA: i32 = 1;
    pub const EINTERP: i32 = 10;
    pub const EBB: i32 = 8;
    pub const EBBPOL: i32 = 9;
    pub const EFBB: i32 = 11;
    pub const ERHO: i32 = 12;
    pub const EPOLSYNCHTH: i32 = 3;
    pub const ESYNCHTHAV: i32 = 4;
    pub const EPOLSYNCHPL: i32 = 2;
    pub const ESYNCHPL: i32 = 6;
    pub const EHYBRIDTH: i32 = 13;
    pub const EHYBRIDPL: i32 = 14;
    pub const EHYBRIDTHPL: i32 = 15;
}

/// Upstream `emis_params`.
#[derive(Clone, Debug, Default)]
pub struct EmisParams {
    pub gmin: Vec<f64>,
    pub mu: Vec<f64>,
    pub nfreq_tab: i32,
    pub freq_tab: Vec<f64>,
}

/// Upstream `emis`.
#[derive(Clone, Debug)]
pub struct Emis {
    pub type_: i32,
    pub neq: usize,
    pub nk: usize,
    pub npts: usize,
    /// emission coefficients (npts x neq): jI, jQ, jU, jV
    pub j: Vec<f64>,
    /// transfer coefficients (npts x nk): alphaI,Q,U,V and Faraday terms
    pub kcoef: Vec<f64>,
    pub rshift: Vec<f64>,
    pub incang: Vec<f64>,
    pub cosne: Vec<f64>,
    pub gmin: Vec<f64>,
    pub tcgs: Vec<f64>,
    pub ncgs: Vec<f64>,
    pub ncgsnth: Vec<f64>,
    pub bcgs: Vec<f64>,
    pub fcol: f64,
    pub p: Vec<f64>,
    pub gmax: f64,
}

impl Emis {
    /// Upstream `select_emissivity_values`.
    pub fn select(ename: &str) -> Emis {
        let mut e = Emis {
            type_: 0,
            neq: 0,
            nk: 0,
            npts: 0,
            j: Vec::new(),
            kcoef: Vec::new(),
            rshift: Vec::new(),
            incang: Vec::new(),
            cosne: Vec::new(),
            gmin: Vec::new(),
            tcgs: Vec::new(),
            ncgs: Vec::new(),
            ncgsnth: Vec::new(),
            bcgs: Vec::new(),
            fcol: 1.8,
            p: Vec::new(),
            gmax: 1e5,
        };
        e.type_ = match ename {
            "POLSYNCHTH" => etype::EPOLSYNCHTH,
            "SYNCHTHAV" => etype::ESYNCHTHAV,
            "POLSYNCHPL" => etype::EPOLSYNCHPL,
            "SYNCHPL" => etype::ESYNCHPL,
            "HYBRIDTH" => etype::EHYBRIDTH,
            "HYBRIDPL" => etype::EHYBRIDPL,
            "HYBRIDTHPL" => etype::EHYBRIDTHPL,
            "lambda" => etype::ELAMBDA,
            "BB" => etype::EBB,
            "FBB" => etype::EFBB,
            "RHO" => etype::ERHO,
            "BBPOL" => etype::EBBPOL,
            "INTERP" => etype::EINTERP,
            "INTERPPOL" => etype::EINTERP,
            _ => {
                eprintln!("WARNING -- Emissivity not recognized: {ename}");
                0
            }
        };
        e.neq = match e.type_ {
            etype::EPOLSYNCHTH
            | etype::EPOLSYNCHPL
            | etype::EBBPOL
            | etype::EHYBRIDTH
            | etype::EHYBRIDPL
            | etype::EHYBRIDTHPL => 4,
            _ => 1,
        };
        e.nk = 1 + e.neq * (e.neq - 1) / 2;
        e
    }

    /// Upstream `initialize_emissivity` (allocation and stored frame data).
    pub fn initialize(&mut self, npts: usize, rshift: &[f64], ang: &[f64], cosne: &[f64]) {
        self.npts = npts;
        self.j = vec![0.0; npts * self.neq];
        self.kcoef = vec![0.0; npts * self.nk];
        self.rshift = rshift.to_vec();
        self.incang = ang.to_vec();
        self.cosne = cosne.to_vec();
        self.gmin = vec![0.0; npts];
        match self.type_ {
            etype::EBB | etype::EFBB | etype::EBBPOL => {
                self.tcgs = vec![0.0; npts];
            }
            etype::EPOLSYNCHPL | etype::ESYNCHPL => {
                self.ncgsnth = vec![0.0; npts];
                self.bcgs = vec![0.0; npts];
                self.p = vec![0.0; npts];
            }
            etype::EPOLSYNCHTH | etype::ESYNCHTHAV => {
                self.ncgs = vec![0.0; npts];
                self.tcgs = vec![0.0; npts];
                self.bcgs = vec![0.0; npts];
            }
            _ => {}
        }
    }

    /// Upstream `assign_emis_params` for the synchrotron types.
    pub fn assign_synch_params(
        &mut self,
        ncgs: &[f64],
        ncgsnth: &[f64],
        bcgs: &[f64],
        tcgs: &[f64],
    ) {
        match self.type_ {
            etype::EPOLSYNCHPL | etype::ESYNCHPL => {
                self.ncgsnth = ncgsnth.to_vec();
                self.bcgs = bcgs.to_vec();
                self.tcgs = tcgs.to_vec();
            }
            etype::EPOLSYNCHTH | etype::ESYNCHTHAV => {
                self.ncgs = ncgs.to_vec();
                self.bcgs = bcgs.to_vec();
                self.tcgs = tcgs.to_vec();
            }
            _ => {}
        }
    }

    /// Upstream `emis_model_synchpl`: store p, gmin, gmax for the power-law
    /// model.
    pub fn set_synchpl_model(&mut self, p: f64, gmin: &[f64], gmax: f64) {
        self.p = vec![p; self.npts];
        self.gmin = gmin.to_vec();
        self.gmax = gmax;
    }

    /// Upstream `assign_emis_params` for the ported types.
    pub fn assign_params(&mut self, tcgs: &[f64]) {
        match self.type_ {
            etype::EBB | etype::EFBB | etype::EBBPOL => {
                self.tcgs = tcgs.to_vec();
            }
            _ => {}
        }
    }

    /// Upstream `calc_emissivity` for the ported types. `nu` is the emitted
    /// frequency at each point (observed frequency divided by the redshift).
    pub fn calc_emissivity(&mut self, nu: &[f64]) {
        let npts = self.npts;
        // K buffer with the upstream 11-column convention:
        // [jI, jQ, jU, jV, alphaI, alphaQ, alphaU, alphaV, rhoQ, rhoU, rhoV]
        let mut kb = vec![0.0f64; npts * 11];
        match self.type_ {
            etype::EBB => {
                for i in 0..npts {
                    kb[i * 11] = bnu(&self.tcgs[i..i + 1], nu[i])[0];
                }
            }
            etype::EFBB => {
                for i in 0..npts {
                    let f = self.fcol;
                    kb[i * 11] = f.powf(-4.0) * bnu(&[(self.tcgs[i] * f)], nu[i])[0];
                }
            }
            etype::EBBPOL => {
                // fbbpolemis: f is hardcoded to 1.8 upstream
                let f = 1.8f64;
                for i in 0..npts {
                    let mut v = f.powf(-4.0) * bnu(&[self.tcgs[i] * f], nu[i])[0];
                    let (interp_i, interp_del) = interp_chandra_tab24_f32(self.cosne[i] as f32);
                    v *= interp_i as f64;
                    kb[i * 11] = v;
                    kb[i * 11 + 1] = v * interp_del as f64;
                }
            }
            etype::ELAMBDA => {
                // lambda(e): e%j=1, e%K=0
                for i in 0..npts {
                    kb[i * 11] = 1.0;
                }
            }
            etype::EPOLSYNCHTH => {
                let out = polsynchth(nu, &self.ncgs, &self.bcgs, &self.tcgs, &self.incang);
                for (i, row) in out.iter().enumerate() {
                    kb[i * 11..i * 11 + 11].copy_from_slice(row);
                }
            }
            etype::ESYNCHTHAV => {
                let out = synchemis(nu, &self.ncgs, &self.bcgs, &self.tcgs);
                for (i, row) in out.iter().enumerate() {
                    kb[i * 11..i * 11 + 11].copy_from_slice(row);
                }
            }
            etype::EPOLSYNCHPL => {
                let out = polsynchpl(
                    nu,
                    &self.ncgsnth,
                    &self.bcgs,
                    &self.incang,
                    &self.p,
                    &self.gmin,
                    self.gmax,
                );
                for (i, row) in out.iter().enumerate() {
                    kb[i * 11..i * 11 + 11].copy_from_slice(row);
                }
            }
            etype::ESYNCHPL => {
                let out = synchpl(
                    nu,
                    &self.ncgsnth,
                    &self.bcgs,
                    &self.incang,
                    &self.p,
                    &self.gmin,
                    self.gmax,
                );
                for (i, row) in out.iter().enumerate() {
                    kb[i * 11..i * 11 + 11].copy_from_slice(row);
                }
            }
            _ => {
                panic!(
                    "calc_emissivity: emissivity type {} not yet ported",
                    self.type_
                );
            }
        }
        // split into j (first neq columns) and K (next nk columns)
        for i in 0..npts {
            for q in 0..self.neq {
                self.j[i * self.neq + q] = kb[i * 11 + q];
            }
            for q in 0..self.nk {
                self.kcoef[i * self.nk + q] = kb[i * 11 + self.neq + q];
            }
        }
    }
}

/// Upstream `rotate_emis`: rotate Q/U emission and transfer coefficients by
/// the comoving-frame angle.
pub fn rotate_emis(e: &mut Emis, s2xi: &[f64], c2xi: &[f64]) {
    let n = e.npts;
    for i in 0..n {
        let jq = e.j[i * e.neq + 1];
        let ju = e.j[i * e.neq + 2];
        let aq = e.kcoef[i * e.nk + 1];
        let au = e.kcoef[i * e.nk + 2];
        let rhoq = e.kcoef[i * e.nk + 4];
        let rhou = e.kcoef[i * e.nk + 5];
        e.j[i * e.neq + 1] = c2xi[i] * jq - s2xi[i] * ju;
        e.j[i * e.neq + 2] = s2xi[i] * jq + c2xi[i] * ju;
        e.kcoef[i * e.nk + 1] = c2xi[i] * aq - s2xi[i] * au;
        e.kcoef[i * e.nk + 2] = s2xi[i] * aq + c2xi[i] * au;
        e.kcoef[i * e.nk + 4] = c2xi[i] * rhoq - s2xi[i] * rhou;
        e.kcoef[i * e.nk + 5] = s2xi[i] * rhoq + c2xi[i] * rhou;
    }
}

/// Upstream `invariant_emis`: scale emission by `g^npow` and absorption by
/// `1/g` (the driver uses `npow=3` for single-point rays and the 2-argument
/// form, `j*=g^2`, `K/=g`, for integrated rays).
pub fn invariant_emis(e: &mut Emis, g: &[f64], npow: i32) {
    let n = e.npts;
    for i in 0..n {
        for q in 0..e.neq {
            e.j[i * e.neq + q] *= g[i].powi(npow);
        }
        for q in 0..e.nk {
            e.kcoef[i * e.nk + q] /= g[i];
        }
    }
}

/// The two-argument upstream `invariant_emis` (`j *= g^2`, `K /= g`).
pub fn invariant_emis_g2(e: &mut Emis, g: &[f64]) {
    let n = e.npts;
    for i in 0..n {
        for q in 0..e.neq {
            e.j[i * e.neq + q] *= g[i] * g[i];
        }
        for q in 0..e.nk {
            e.kcoef[i * e.nk + q] /= g[i];
        }
    }
}

/// Upstream `assign_emis_params` (free function form matching the driver).
pub fn assign_emis_params(e: &mut Emis, tcgs: &[f64]) {
    e.assign_params(tcgs);
}

/// Upstream `calc_emissivity`.
pub fn calc_emissivity(nu: &[f64], e: &mut Emis, _ep: &EmisParams) {
    e.calc_emissivity(nu);
}
