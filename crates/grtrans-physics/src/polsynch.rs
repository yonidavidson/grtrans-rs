//! Polarized synchrotron emissivities.
//!
//! Translation target: `polsynchemis.f90` (upstream GRTRANS).
//! Currently ported: the Planck function (`bnu`), the power-law
//! coefficients (`polsynchpl`, `synchpl`) and their interpolation tables
//! (extracted verbatim into [`crate::polsynch_tables`]).

use crate::polsynch_tables::TABLES;
use grtrans_core::constants::{C, C2, E as EC, H, K, M as ME, PI};
use std::sync::OnceLock;

/// Planck spectrum for an array of temperatures (cgs), upstream `bnu`
/// (`polsynchemis.f90` lines 1014-1032).
///
/// Uses the Rayleigh-Jeans limit for `h nu / k T < 1e-6` and replaces exact
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

/// The log-transformed interpolation tables (upstream applies `log` at
/// initialization).
pub struct Tables {
    pub gxvals: [f64; 181],
    pub gyvals: [f64; 543],
    pub xvals: [f64; 81],
    pub yavals: [f64; 243],
    pub ypvals: [f64; 243],
    pub yvvals: [f64; 243],
    pub yapvals: [f64; 243],
    pub yavvals: [f64; 243],
}

pub fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let mut t = Tables {
            gxvals: [0.0; 181],
            gyvals: [0.0; 543],
            xvals: [0.0; 81],
            yavals: [0.0; 243],
            ypvals: [0.0; 243],
            yvvals: [0.0; 243],
            yapvals: [0.0; 243],
            yavvals: [0.0; 243],
        };
        for (dst, src) in t.gxvals.iter_mut().zip(TABLES.gxvals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.gyvals.iter_mut().zip(TABLES.gyvals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.xvals.iter_mut().zip(TABLES.xvals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.yavals.iter_mut().zip(TABLES.yavals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.ypvals.iter_mut().zip(TABLES.ypvals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.yvvals.iter_mut().zip(TABLES.yvvals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.yapvals.iter_mut().zip(TABLES.yapvals.iter()) {
            *dst = src.ln();
        }
        for (dst, src) in t.yavvals.iter_mut().zip(TABLES.yavvals.iter()) {
            *dst = src.ln();
        }
        t
    })
}

/// Upstream `find_inds_gxp`: table indices for one point.
pub fn find_inds_gxp(x: f64, p: f64, xv: &[f64]) -> (usize, usize) {
    let nx = xv.len();
    let xx = x.ln();
    let mut ix = ((xx - xv[0]) / (xv[nx - 1] - xv[0]) * nx as f64) as i64;
    // merge(merge(ix,1,ix.gt.1),nx-1,ix.lt.(nx-1))
    if ix <= 1 {
        ix = 1;
    }
    if ix >= (nx as i64 - 1) {
        ix = nx as i64 - 1;
    }
    let ix = ix as usize; // 1-based
                          // iy: p<=3.2 -> 0; 3.2<p<5 -> 1; p>=5 -> 2
    let iy = if p <= 3.2 {
        0
    } else if p < 5.0 {
        1
    } else {
        2
    };
    (ix, iy)
}

/// Upstream `interp_gxp`: interpolate in log space and exponentiate.
pub fn interp_gxp(x: f64, ix: usize, iy: usize, xv: &[f64], yv: &[f64]) -> f64 {
    let xx = x.ln();
    let indx = iy * xv.len() + (ix - 1); // 0-based
    let yix = yv[indx];
    let yix1 = yv[indx + 1];
    let xix = xv[ix - 1];
    let xix1 = xv[ix];
    let slope = (yix1 - yix) / (xix1 - xix);
    (yix + slope * (xx - xix)).exp()
}

/// The six G-functions of `get_polsynchpl_facs`.
pub fn get_polsynchpl_facs(
    xmin: &[f64],
    xmax: &[f64],
    p: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let t = tables();
    let n = xmin.len();
    let mut gfac = vec![0.0; n];
    let mut gpfac = vec![0.0; n];
    let mut gvfac = vec![0.0; n];
    let mut gafac = vec![0.0; n];
    let mut gapfac = vec![0.0; n];
    let mut gavfac = vec![0.0; n];
    for i in 0..n {
        let (ix2, iy2) = find_inds_gxp(xmax[i], p[i], &t.gxvals);
        let (ix1, iy1) = find_inds_gxp(xmin[i], p[i], &t.gxvals);
        gfac[i] = interp_gxp(xmax[i], ix2, iy2, &t.gxvals, &t.gyvals)
            - interp_gxp(xmin[i], ix1, iy1, &t.gxvals, &t.gyvals);
        let (ix2, iy2) = find_inds_gxp(xmax[i], p[i], &t.xvals);
        let (ix1, iy1) = find_inds_gxp(xmin[i], p[i], &t.xvals);
        gafac[i] = interp_gxp(xmax[i], ix2, iy2, &t.xvals, &t.yavals)
            - interp_gxp(xmin[i], ix1, iy1, &t.xvals, &t.yavals);
        gpfac[i] = interp_gxp(xmax[i], ix2, iy2, &t.xvals, &t.ypvals)
            - interp_gxp(xmin[i], ix1, iy1, &t.xvals, &t.ypvals);
        gapfac[i] = interp_gxp(xmax[i], ix2, iy2, &t.xvals, &t.yapvals)
            - interp_gxp(xmin[i], ix1, iy1, &t.xvals, &t.yapvals);
        gvfac[i] = interp_gxp(xmax[i], ix2, iy2, &t.xvals, &t.yvvals)
            - interp_gxp(xmin[i], ix1, iy1, &t.xvals, &t.yvvals);
        gavfac[i] = interp_gxp(xmax[i], ix2, iy2, &t.xvals, &t.yavvals)
            - interp_gxp(xmin[i], ix1, iy1, &t.xvals, &t.yavvals);
    }
    (gfac, gpfac, gvfac, gafac, gapfac, gavfac)
}

/// Upstream `polsynchpl`: polarized power-law synchrotron coefficients
/// (11 outputs per point).
#[allow(clippy::too_many_arguments)]
pub fn polsynchpl(
    nu: &[f64],
    nnth: &[f64],
    b: &[f64],
    th: &[f64],
    p: &[f64],
    gmin: &[f64],
    gmax: f64,
) -> Vec<[f64; 11]> {
    let n = nu.len();
    let mut out = vec![[0.0f64; 11]; n];
    let nufloor = 1e-10;
    let thsafe = 1e-10;
    let mut xmin = vec![0.0; n];
    let mut xmax = vec![0.0; n];
    let mut tanth = vec![0.0; n];
    let mut sinth = vec![0.0; n];
    let mut nubperp = vec![0.0; n];
    let mut a_pref = vec![0.0; n];
    for i in 0..n {
        tanth[i] = th[i].tan() + if th[i].cos() >= 0.0 { 1.0 } else { -1.0 } * thsafe;
        sinth[i] = th[i].sin() + thsafe;
        nubperp[i] = EC * b[i] / ME / C / 2.0 / PI * sinth[i] + nufloor;
        let nucmin = 1.5 * nubperp[i] * gmin[i].powi(2);
        let nucmax = 1.5 * nubperp[i] * gmax * gmax;
        xmin[i] = nu[i] / nucmin;
        xmax[i] = nu[i] / nucmax;
        a_pref[i] = (p[i] - 1.0) * nnth[i] / (gmin[i].powf(1.0 - p[i]) - gmax.powf(1.0 - p[i]));
    }
    let (gxfac, gpfac, gvfac, gafac, gapfac, gavfac) = get_polsynchpl_facs(&xmin, &xmax, p);
    for i in 0..n {
        let omega0 = nubperp[i] * 2.0 * PI;
        let omega = nu[i] * 2.0 * PI;
        let a = a_pref[i];
        let jfac = a * EC * EC / C * 3.0f64.sqrt() / 4.0
            * (3.0 * nubperp[i] / 2.0 / nu[i]).powf((p[i] - 1.0) / 2.0)
            * nubperp[i];
        let ji = jfac * gxfac[i];
        let jq = jfac * gpfac[i];
        let jv = jfac * 4.0 / 3.0 / tanth[i] * (3.0 * nubperp[i] / 2.0 / nu[i]).sqrt() * gvfac[i];
        let alpha = (p[i] - 1.0) / 2.0;
        let kperp = a * EC * EC / ME / C / nubperp[i];
        let nui = gmin[i] * gmin[i] * nubperp[i];
        let kstaralphaq = 1.0;
        let kstaralphav = 2.0 * (alpha + 3.0 / 2.0) / (alpha + 1.0);
        let kstarq = kstaralphaq
            * kperp
            * (nubperp[i] / nu[i]).powi(3)
            * gmin[i].powf(-2.0 * alpha + 1.0)
            * (1.0 - (nui / nu[i]).powf(alpha - 1.0 / 2.0))
            * (alpha - 1.0 / 2.0).powi(-1);
        let kstarv = kstaralphav
            * kperp
            * (nubperp[i] / nu[i]).powi(2)
            * gmin[i].ln()
            * gmin[i].powf(-2.0 * (alpha + 1.0))
            / tanth[i];
        let afac = (2.0 * PI).powi(3) * a * EC * EC * 3.0f64.sqrt() * omega0 * (p[i] + 2.0)
            / 32.0
            / PI.powi(2)
            / ME
            / C
            / omega.powi(2)
            * (2.0 * omega / 3.0 / omega0).powf(-p[i] / 2.0);
        let ai = afac * gafac[i];
        let aq = afac * gapfac[i];
        let av =
            afac * 4.0 / 3.0 / tanth[i] * gavfac[i] * (2.0 * omega / 3.0 / omega0).powf(-1.0 / 2.0);
        out[i][0] = ji;
        out[i][1] = jq;
        out[i][2] = 0.0;
        out[i][3] = jv;
        out[i][4] = ai;
        out[i][5] = aq;
        out[i][6] = 0.0;
        out[i][7] = av;
        out[i][8] = kstarq;
        out[i][9] = 0.0;
        out[i][10] = kstarv;
    }
    out
}

/// Upstream `synchpl`: unpolarized power-law synchrotron coefficients.
#[allow(clippy::too_many_arguments)]
pub fn synchpl(
    nu: &[f64],
    nnth: &[f64],
    b: &[f64],
    th: &[f64],
    p: &[f64],
    gmin: &[f64],
    gmax: f64,
) -> Vec<[f64; 11]> {
    let n = nu.len();
    let mut out = vec![[0.0f64; 11]; n];
    let t = tables();
    let nufloor = 1e-10;
    for i in 0..n {
        let sinth = th[i].sin() + nufloor;
        let nubperp = EC * b[i] / ME / C / 2.0 / PI * sinth + nufloor;
        let nucmin = 1.5 * nubperp * gmin[i].powi(2);
        let nucmax = 1.5 * nubperp * gmax * gmax;
        let omega0 = nubperp * 2.0 * PI;
        let omega = nu[i] * 2.0 * PI;
        let xmin = nu[i] / nucmin;
        let xmax = nu[i] / nucmax;
        let a = (p[i] - 1.0) * nnth[i] / (gmin[i].powf(1.0 - p[i]) - gmax.powf(1.0 - p[i]));
        let (ix2, iy2) = find_inds_gxp(xmax, p[i], &t.gxvals);
        let (ix1, iy1) = find_inds_gxp(xmin, p[i], &t.gxvals);
        let gfac = interp_gxp(xmax, ix2, iy2, &t.gxvals, &t.gyvals)
            - interp_gxp(xmin, ix1, iy1, &t.gxvals, &t.gyvals);
        let (ix1, iy1) = find_inds_gxp(xmin, p[i], &t.xvals);
        let (ix2, iy2) = find_inds_gxp(xmax, p[i], &t.xvals);
        let gafac = interp_gxp(xmax, ix2, iy2, &t.xvals, &t.yavals)
            - interp_gxp(xmin, ix1, iy1, &t.xvals, &t.yavals);
        let jfac = a * EC * EC / C * 3.0f64.sqrt() / 4.0
            * (3.0 * nubperp / 2.0 / nu[i]).powf((p[i] - 1.0) / 2.0)
            * nubperp;
        let ji = jfac * gfac;
        let afac = (2.0 * PI).powi(3) * a * EC * EC * 3.0f64.sqrt() * omega0 * (p[i] + 2.0)
            / 32.0
            / PI.powi(2)
            / ME
            / C
            / omega.powi(2)
            * (2.0 * omega / 3.0 / omega0).powf(-p[i] / 2.0);
        let ai = afac * gafac;
        out[i][0] = ji;
        out[i][4] = ai;
    }
    out
}

// ---------------------------------------------------------------------------
// Thermal synchrotron (upstream polsynchth, synchemis and their fitting
// functions; Huang et al. 2009, Shcherbakov 2008, Mahadevan et al. 1996).
// ---------------------------------------------------------------------------

/// Upstream `shffunc` (Shcherbakov 2008 fitting function F(X)).
pub fn shffunc(x: f64) -> f64 {
    2.011 * (-x.powf(1.035) / 4.7).exp()
        - (x / 2.0).cos() * (-x.powf(1.2) / 2.73).exp()
        - 0.011 * (-x / 47.2).exp()
}

/// Upstream `jffunc` (modified F(X) matching the Jones & Hardee limit).
pub fn jffunc(x: f64) -> f64 {
    let extraterm = (0.011 * (-x / 47.2).exp()
        - 2f64.powf(-1.0 / 3.0) / 3f64.powf(23.0 / 6.0) * PI * 1e4 * (x + 1e-16).powf(-8.0 / 3.0))
        * (0.5 + 0.5 * ((x.ln() - 120f64.ln()) / 0.1).tanh());
    2.011 * (-x.powf(1.035) / 4.7).exp()
        - (x / 2.0).cos() * (-x.powf(1.2) / 2.73).exp()
        - 0.011 * (-x / 47.2).exp()
        + extraterm
}

/// Upstream `shgmfunc` (modified Shcherbakov 2008 G(X)).
pub fn shgmfunc(x: f64) -> f64 {
    0.43793091 * (1.0 + 0.00185777 * x.powf(1.50316886)).ln()
}

/// Upstream `iix` (Mahadevan et al. 1996 fitting function).
pub fn iix(x: f64) -> f64 {
    2.5651
        * (1.0 + 1.92 / x.powf(1.0 / 3.0) + 0.9977 / x.powf(2.0 / 3.0))
        * (-1.8899 * x.powf(1.0 / 3.0)).exp()
}

/// Upstream `iqx`.
pub fn iqx(x: f64) -> f64 {
    2.5651
        * (1.0 + 0.93193 / x.powf(1.0 / 3.0) + 0.499873 / x.powf(2.0 / 3.0))
        * (-1.8899 * x.powf(1.0 / 3.0)).exp()
}

/// Upstream `ivx`.
pub fn ivx(x: f64) -> f64 {
    (1.81384 / x + 3.42319 / x.powf(2.0 / 3.0) + 0.0292545 / x.sqrt() + 2.03773 / x.powf(1.0 / 3.0))
        * (-1.8899 * x.powf(1.0 / 3.0)).exp()
}

/// Upstream `ipx` (Mahadevan et al. 1998 fitting function).
pub fn ipx(x: f64) -> f64 {
    4.0505 / (x.ln() / 6.0).exp()
        * (1.0 + 0.40 / x.sqrt().sqrt() + 0.5316 / x.sqrt())
        * (-1.8899 * (x.ln() / 3.0).exp()).exp()
}

/// Upstream `polsynchth`: polarized thermal synchrotron coefficients
/// (11 outputs per point).
pub fn polsynchth(nu: &[f64], n: &[f64], b: &[f64], t: &[f64], theta: &[f64]) -> Vec<[f64; 11]> {
    use grtrans_core::bessel::{besselk, besselk0, besselk1};
    let np = nu.len();
    let mut out = vec![[0.0f64; 11]; np];
    let thetaemin = 1e-10;
    for i in 0..np {
        let thetae = K * t[i] / ME / C / C + thetaemin;
        let sin_th = theta[i].sin();
        let nuc = 3.0 * EC * b[i] * sin_th / 4.0 / PI / ME / C * thetae * thetae + 1.0;
        let xm = nu[i] / nuc;
        let ji = EC * EC / C / 3.0f64.sqrt() / 2.0 * n[i] / (thetae * thetae) * nu[i] * iix(xm);
        let jq = EC * EC / C / 3.0f64.sqrt() / 2.0 * n[i] / (thetae * thetae) * nu[i] * iqx(xm);
        let jv =
            4.0 * EC * EC / C / 3.0 / 3.0f64.sqrt() / theta[i].tan() * n[i] / 2.0 / thetae.powi(3)
                * nu[i]
                * ivx(xm);
        let ju = 0.0;
        let bnutnu = bnu(&[t[i]], nu[i])[0];
        let ai = ji / bnutnu;
        let aq = jq / bnutnu;
        let av = jv / bnutnu;
        let au = ju / bnutnu;
        let rhou = au;
        let wp2 = 4.0 * PI * n[i] * EC * EC / ME;
        let omega0 = EC * b[i] / ME / C;
        let xarg = thetae * (2.0f64.sqrt() * sin_th * (1e3 * omega0 / 2.0 / PI / nu[i])).sqrt();
        let (eps11m22, eps12);
        if thetae > 1e-2 {
            let k1k2 = besselk1(1.0 / thetae) / besselk(2, 1.0 / thetae);
            eps11m22 = jffunc(xarg) * wp2 * omega0 * omega0 / (2.0 * PI * nu[i]).powi(4)
                * (k1k2 + 6.0 * thetae)
                * sin_th
                * sin_th;
            let step = 0.5 + 0.5 * ((thetae - 1.0) / 0.05).tanh();
            eps12 = wp2 * omega0 / (2.0 * PI * nu[i]).powi(3)
                * (besselk0(1.0 / thetae) - step * shgmfunc(xarg))
                / besselk(2, 1.0 / thetae)
                * theta[i].cos();
        } else {
            eps11m22 = jffunc(xarg) * wp2 * omega0 * omega0 / (2.0 * PI * nu[i]).powi(4)
                * (1.0 + 6.0 * thetae)
                * sin_th
                * sin_th;
            eps12 = wp2 * omega0 / (2.0 * PI * nu[i]).powi(3) * theta[i].cos();
        }
        let rhov = 2.0 * PI * nu[i] / C * eps12;
        let rhoq = 2.0 * PI * nu[i] / 2.0 / C * eps11m22;
        out[i] = [ji, jq, ju, jv, ai, aq, au, av, rhoq, rhou, rhov];
    }
    out
}

/// Upstream `synchemis` (Mahadevan et al. 1996 unpolarized thermal
/// synchrotron; used by `SYNCHTHAV`).
pub fn synchemis(nu: &[f64], n: &[f64], b: &[f64], t: &[f64]) -> Vec<[f64; 11]> {
    let np = nu.len();
    let mut out = vec![[0.0f64; 11]; np];
    let thetaemin = 1e-10;
    for i in 0..np {
        let thetae = K * t[i] / ME / C2 + thetaemin;
        let nucrit = 3.0 * EC * b[i] / (4.0 * PI * ME * C) * thetae * thetae + 1.0;
        let xm = nu[i] / nucrit;
        let jnu = 4.43e-30 / 2.0 * nu[i] * n[i] * ipx(xm) / (thetae * thetae);
        let anu = if jnu.abs() > 0.0 {
            jnu / bnu(&[t[i]], nu[i])[0]
        } else {
            0.0
        };
        out[i][0] = jnu;
        out[i][4] = anu;
    }
    out
}
