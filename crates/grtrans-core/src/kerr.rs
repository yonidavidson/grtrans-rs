//! Kerr spacetime geometry, frames, geodesic constants and polarization
//! transport.
//!
//! Direct translation of `kerr.f90` (upstream GRTRANS). Functions operate in
//! geometrized units (G = c = M = 1) with Boyer-Lindquist (BL) coordinates
//! `x^mu = (t, r, theta, phi)` unless stated otherwise.
//!
//! Metric component ordering follows upstream: the ten components are the
//! upper triangle `[g00, g01, g02, g03, g11, g12, g13, g22, g23, g33]`.

use crate::constants::{G, MP, MSUN, PI, SIGT};
use crate::four_vector::{FourVector, Metric};
use std::f64::consts::PI as PI_F64;

/// Boyer-Lindquist covariant metric components (upstream `blmetric_cov`).
pub fn blmetric_cov(r: f64, th: f64, a: f64) -> Metric {
    let cth = th.cos();
    let sth = th.sin();
    let delta = r * r - 2.0 * r + a * a;
    let rho2 = r * r + a * a * cth * cth;
    let sigma = (r * r + a * a).powi(2) - a * a * delta * sth * sth;
    let mut m = [0.0f64; 10];
    m[0] = -(delta - a * a * sth * sth) / rho2;
    m[3] = -2.0 * a * r * sth * sth / rho2;
    m[4] = rho2 / delta;
    m[7] = rho2;
    m[9] = sigma / rho2 * sth * sth;
    Metric(m)
}

/// Boyer-Lindquist contravariant metric components (upstream `blmetric_con`).
pub fn blmetric_con(r: f64, th: f64, a: f64) -> Metric {
    let cth = th.cos();
    let sth = th.sin();
    let delta = r * r - 2.0 * r + a * a;
    let rho2 = r * r + a * a * cth * cth;
    let mut m = [0.0f64; 10];
    m[0] = -((r * r + a * a).powi(2) - a * a * delta * sth * sth) / rho2 / delta;
    m[3] = -2.0 * a * r / rho2 / delta;
    m[4] = delta / rho2;
    m[7] = 1.0 / rho2;
    m[9] = (delta - a * a * sth * sth) / delta / rho2 / sth / sth;
    Metric(m)
}

/// Boyer-Lindquist covariant metric computed in single precision and
/// widened to f64 (upstream `blmetric_cov_real`). Used by the thin-disk
/// fluid model, which passes `real` arguments and receives an `real` array.
pub fn blmetric_cov_f32(r: f32, th: f32, a: f32) -> Metric {
    let cth = th.cos();
    let sth = th.sin();
    let delta = r * r - 2.0 * r + a * a;
    let rho2 = r * r + a * a * cth * cth;
    let sigma = (r * r + a * a).powi(2) - a * a * delta * sth * sth;
    let mut m = [0.0f64; 10];
    m[0] = (-(delta - a * a * sth * sth) / rho2) as f64;
    m[3] = (-2.0 * a * r * sth * sth / rho2) as f64;
    m[4] = (rho2 / delta) as f64;
    m[7] = rho2 as f64;
    m[9] = (sigma / rho2 * sth * sth) as f64;
    Metric(m)
}

/// Kerr-Schild spherical covariant metric (upstream `ksmetric_cov`).
pub fn ksmetric_cov(r: f64, th: f64, _ph: f64, a: f64) -> Metric {
    let ctheta = th.cos();
    let stheta = th.sin();
    let rho2 = r * r + a * a * ctheta * ctheta;
    let psi4 = 2.0 * r / rho2;
    let mut m = [0.0f64; 10];
    m[0] = -(1.0 - psi4);
    m[1] = psi4;
    m[3] = -a * stheta * stheta * psi4;
    m[4] = 1.0 + psi4;
    m[6] = -a * stheta * stheta * (1.0 + psi4);
    m[7] = rho2;
    m[9] = stheta * stheta * (rho2 + a * a * (1.0 + 2.0 * r / rho2) * stheta * stheta);
    Metric(m)
}

/// Roche/ISCO radius from Bardeen, Press & Teukolsky, upstream `calc_rms`.
pub fn calc_rms(a: f64) -> f64 {
    let a = a as f32; // upstream uses default `real` (f32) here
    let z1 = 1.0
        + (1.0 - a * a).powf(1.0 / 3.0) * ((1.0 + a).powf(1.0 / 3.0) + (1.0 - a).powf(1.0 / 3.0));
    let z2 = (3.0 * a * a + z1 * z1).sqrt();
    let s = ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt();
    // Fortran sign(s, a) = |s| * sign(a), with sign(+0.0) = +1
    let signed = if a >= 0.0 { s } else { -s };
    (3.0 + z2 - signed) as f64
}

/// Energy, angular momentum and radius of the ISCO (upstream
/// `calc_rms_constants`).
pub fn calc_rms_constants(a: f64) -> (f64, f64, f64) {
    let rms = calc_rms(a);
    let v = 1.0 / rms.sqrt();
    let ems =
        (1.0 - 2.0 * v * v + a * v * v * v) / (1.0 - 3.0 * v * v + 2.0 * a * v * v * v).sqrt();
    let lms = rms * v * (1.0 - 2.0 * a * v * v * v + a * a * v * v * v * v)
        / (1.0 - 3.0 * v * v + 2.0 * a * v * v * v).sqrt();
    (ems, lms, rms)
}

/// Krolik & Hawley (2002) correction factor for the thin-disk flux
/// (upstream `krolikc`), evaluated at a single radius.
pub fn krolikc(r: f64, a: f64) -> f64 {
    let a = a as f32;
    let r = r as f32;
    let pi = std::f32::consts::PI;
    let rms = calc_rms(a as f64) as f32;
    let y = r.sqrt();
    let yms = rms.sqrt();
    let y1 = 2.0 * ((a.acos() - pi) / 3.0).cos();
    let y2 = 2.0 * ((a.acos() + pi) / 3.0).cos();
    let y3 = -2.0 * (a.acos() / 3.0).cos();
    let arg1 = 3.0 * a / (2.0 * y);
    let arg2 = 3.0 * (y1 - a).powi(2) / (y * y1 * (y1 - y2) * (y1 - y3));
    let arg3 = 3.0 * (y2 - a).powi(2) / (y * y2 * (y2 - y1) * (y2 - y3));
    let arg4 = 3.0 * (y3 - a).powi(2) / (y * y3 * (y3 - y1) * (y3 - y2));
    (1.0 - yms / y
        - arg1 * (y / yms).ln()
        - arg2 * ((y - y1) / (yms - y1)).ln()
        - arg3 * ((y - y2) / (yms - y2)).ln()
        - arg4 * ((y - y3) / (yms - y3)).ln()) as f64
}

/// Eddington luminosity for mass `m` in solar masses (upstream `ledd`).
pub fn ledd(m: f32) -> f64 {
    4.0 * PI * G * (m as f64) * MSUN * MP * crate::constants::C / SIGT
}

/// Convert t or phi from Boyer-Lindquist to Kerr-Schild (upstream
/// `bl2ks_time`).
pub fn bl2ks_time(r: f64, x: f64, a: f64) -> f64 {
    x + (r * r - 2.0 * r + a * a).ln()
        + 1.0 / (2.0 * (1.0 - a * a).sqrt())
            * ((r - 1.0 - (1.0 - a * a).sqrt()) / (r - 1.0 + (1.0 - a * a).sqrt())).ln()
}

/// Convert phi from Boyer-Lindquist to Kerr-Schild (upstream `bl2ks_phi`).
pub fn bl2ks_phi(r: f64, x: f64, a: f64) -> f64 {
    x + a / (2.0 * (1.0 - a * a).sqrt())
        * ((r - 1.0 - (1.0 - a * a).sqrt()) / (r - 1.0 + (1.0 - a * a).sqrt())).ln()
}

/// Convert a Kerr-Schild spherical 4-velocity to Boyer-Lindquist
/// (upstream `uks2ubl`, from Font et al. 1999). `r` is the BL radius and
/// `a` the spin.
pub fn uks2ubl(fut: &FourVector, r: f64, a: f64) -> FourVector {
    let delta = r * r - 2.0 * r + a * a;
    let mut fu = *fut;
    fu.data[0] = fut.data[0] - 2.0 * r / delta * fut.data[1];
    fu.data[3] = fut.data[3] - a / delta * fut.data[1];
    fu
}

/// LNRF components of a velocity, upstream `lnrf_frame`.
pub fn lnrf_frame(vr: f64, vt: f64, omega: f64, r: f64, a: f64, th: f64) -> (f64, f64, f64) {
    let mu = th.cos();
    let d = r * r - 2.0 * r + a * a;
    let ar = (r * r + a * a).powi(2) - a * a * d * (1.0 - mu * mu);
    let rho = r * r + a * a * mu * mu;
    let enu = (d * rho / ar).sqrt();
    let emu1 = (rho / d).sqrt();
    let emu2 = rho.sqrt();
    let epsi = (1.0 - mu * mu).sqrt() * (ar / rho).sqrt();
    let om = 2.0 * a * r / ar;
    let mut vrl = emu1 / enu * vr;
    let mut vtl = emu2 / enu * vt;
    let mut vpl = epsi / enu * (omega - om);
    if !(d > 0.0) {
        vrl = 0.0;
        vtl = 0.0;
        vpl = 0.0;
    }
    (vrl, vtl, vpl)
}

/// Inverse LNRF transform, upstream `lnrf_frame_inv`.
pub fn lnrf_frame_inv(vr: f64, vt: f64, omega: f64, r: f64, a: f64, th: f64) -> (f64, f64, f64) {
    let mu = th.cos();
    let d = r * r - 2.0 * r + a * a;
    let ar = (r * r + a * a).powi(2) - a * a * d * (1.0 - mu * mu);
    let rho = r * r + a * a * mu * mu;
    let enu = (d * rho / ar).sqrt();
    let emu1 = (rho / d).sqrt();
    let emu2 = rho.sqrt();
    let epsi = (1.0 - mu * mu).sqrt() * (ar / rho).sqrt();
    let om = 2.0 * a * r / ar;
    let mut vrl = enu / emu1 * vr;
    let mut vtl = enu / emu2 * vt;
    let mut vpl = enu / epsi * omega + om;
    if !(d > 0.0) {
        vrl = 0.0;
        vtl = 0.0;
        vpl = 0.0;
    }
    (vrl, vtl, vpl)
}

/// Inverse LNRF transform in single precision (upstream
/// `lnrf_frame_inv_real`; selected by gfortran for FFJET's mixed-kind call).
pub fn lnrf_frame_inv_f32(
    vr: f32,
    vt: f32,
    omega: f32,
    r: f32,
    a: f32,
    th: f32,
) -> (f32, f32, f32) {
    let mu = th.cos();
    let d = r * r - 2.0 * r + a * a;
    let ar = (r * r + a * a).powi(2) - a * a * d * (1.0 - mu * mu);
    let rho = r * r + a * a * mu * mu;
    let enu = (d * rho / ar).sqrt();
    let emu1 = (rho / d).sqrt();
    let emu2 = rho.sqrt();
    let epsi = (1.0 - mu * mu).sqrt() * (ar / rho).sqrt();
    let om = 2.0 * a * r / ar;
    let mut vrl = enu / emu1 * vr;
    let mut vtl = enu / emu2 * vt;
    let mut vpl = enu / epsi * omega + om;
    if !(d > 0.0) {
        vrl = 0.0;
        vtl = 0.0;
        vpl = 0.0;
    }
    (vrl, vtl, vpl)
}

/// Covariant wave vector from the constants of motion (upstream
/// `calc_nullp`).
///
/// * `q2`: Carter constant; `l`: axial angular momentum; `a`: spin.
/// * `r`, `mu = cos(theta)`: position.
/// * `su`, `smu`: signs of the radial and polar motion.
/// * `rcomp`/`thcomp`: when true selects the alternative branch used by the
///   geokerr ray reconstruction (both are passed as present upstream).
///
/// Negative radicands are clamped to zero exactly as upstream `merge`.
pub fn calc_nullp(
    q2: f64,
    l: f64,
    a: f64,
    r: f64,
    mu: f64,
    su: f64,
    smu: f64,
    rcomp: bool,
    thcomp: bool,
) -> [f64; 4] {
    let u = 1.0 / r;
    let rho2 = r * r + a * a * mu * mu;
    let delta = r * r - 2.0 * r + a * a;

    let mut mfunc = q2 + (a * a - q2 - l * l) * mu * mu - a * a * mu * mu * mu * mu;
    if !(mfunc > 0.0) {
        mfunc = 0.0;
    }
    let pmu = if thcomp {
        // negative of RB94 to account for dlambda -> -dlambda for forward
        // in time (upstream note 12/11/12)
        1.0 / rho2 * smu * (mfunc / (1.0 - mu * mu)).sqrt()
    } else {
        -smu * mfunc.sqrt() / rho2
    };

    let mut ufunc = 1.0 + (a * a - q2 - l * l) * u * u + 2.0 * ((a - l) * (a - l) + q2) * u * u * u
        - a * a * q2 * u * u * u * u;
    if !(ufunc > 0.0) {
        ufunc = 0.0;
    }
    let pu = if rcomp {
        let rfunc = r * r * ufunc.sqrt();
        su * rfunc / rho2
    } else {
        -su * ufunc.sqrt() / rho2
    };

    let pt =
        (-a * (a * (1.0 - mu * mu) - l) + (r * r + a * a) / delta * (r * r + a * a - a * l)) / rho2;
    let pphi = (-a + l / (1.0 - mu * mu) + a / delta * (r * r + a * a - a * l)) / rho2;
    [pt, pu, pmu, pphi]
}

/// Redshift factor for a Keplerian thin disk (upstream `calcg`), evaluated
/// at one point. `tpm`, `tpr` are the theta/phi turning-point parities.
#[allow(clippy::too_many_arguments)]
pub fn calcg(
    u: f64,
    mu: f64,
    q2: f64,
    l: f64,
    a: f64,
    tpm: i32,
    tpr: i32,
    su: f64,
    sm: f64,
    vrl: f64,
    vtl: f64,
    vpl: f64,
) -> f64 {
    let one = 1.0f64;
    let two = 2.0f64;
    let r = one / u;
    let z1 = one
        + (one - a * a).powf(1.0 / 3.0) * ((one + a).powf(1.0 / 3.0) + (one - a).powf(1.0 / 3.0));
    let z2 = (3.0 * a * a + z1 * z1).sqrt();
    // upstream computes RMS here (unused); kept for transliteration fidelity
    let inner = ((3.0 - z1) * (3.0 + z1 + two * z2)).sqrt();
    let signed_inner = if a >= 0.0 { inner } else { -inner };
    let _rms = 3.0 + z2 - signed_inner;
    let d = r * r - two * r + a * a;
    let ar = (r * r + a * a).powi(2) - a * a * d * (one - mu * mu);
    let rho = r * r + a * a * mu * mu;
    let enu = (d * rho / ar).sqrt();
    let emu1 = (rho / d).sqrt();
    let emu2 = rho.sqrt();
    let epsi = (one - mu * mu).sqrt() * (ar / rho).sqrt();
    let om = two * a * r / ar;
    // (-1)**tpr and (-1)**tpm with Fortran integer-exponent semantics
    let sr = if tpr % 2 == 0 { su } else { -su };
    let st = if tpm % 2 == 0 { -sm } else { sm };
    let mut omega = enu / epsi * vpl + om;
    if epsi == 0.0 {
        omega = 0.0;
    }
    let gam = one / (one - (vrl * vrl + vtl * vtl + vpl * vpl)).sqrt();
    let mut rr = -a * a * q2 * u.powi(4)
        + two * u.powi(3) * (q2 + (a - l).powi(2))
        + u * u * (a * a - q2 - l * l)
        + one;
    if !(rr >= 0.0) {
        rr = 0.0;
    }
    let rr = rr.sqrt() * r * r;
    let mut tt = (q2 + mu * mu * (a * a - l * l - q2) - a * a * mu.powi(4)) / (one - mu * mu);
    if !(tt >= 0.0) {
        tt = 0.0;
    }
    let tt = tt.sqrt();
    enu / gam
        / (one - l * omega - emu1 * enu * vrl / rho * sr * rr - emu2 * enu * vtl / rho * st * tt)
}

/// Angle between the wave vector and magnetic field in the fluid frame
/// (upstream `calc_kb_ang`, Broderick 2004 covariant method).
pub fn calc_kb_ang(k: &FourVector, b: &FourVector, u: &FourVector, r: f64, th: f64, a: f64) -> f64 {
    let tmetric = blmetric_cov(r, th, a); // symmetric: transpose == self
    let mut b = *b;
    let mut u = *u;
    let mut k = *k;
    b.assign_metric(tmetric);
    u.assign_metric(tmetric);
    k.assign_metric(tmetric);
    let bdotk = b.dot(&k);
    let bdotb = b.dot(&b);
    let kdotk = k.dot(&k);
    let om = -k.dot(&u);
    let mut cdot2 = bdotk * bdotk / bdotb / (kdotk + om * om);
    if cdot2 < 0.0 {
        cdot2 = 0.0;
    }
    if cdot2 > 1.0 {
        cdot2 = 1.0;
    }
    let sign = if bdotk >= 0.0 { 1.0 } else { -1.0 };
    (sign * cdot2.sqrt()).acos()
}

/// Walker-Penrose polarization transport (upstream `transport_perpk`).
///
/// Returns `(f1, f2, f3)`, the components of the parallel-transported
/// polarization vector in BL coordinates assuming `f0 = 0`.
pub fn transport_perpk(
    kvec: &FourVector,
    r: f64,
    th: f64,
    a: f64,
    metric: &Metric,
    kap1: f64,
    kap2: f64,
) -> (f64, f64, f64) {
    let g03 = metric.0[3];
    let g11 = metric.0[4];
    let g22 = metric.0[7];
    let g33 = metric.0[9];
    let cth = th.cos();
    let sth = th.sin();
    let k0 = kvec.data[0];
    let k1 = kvec.data[1];
    let k2 = kvec.data[2];
    let k3 = kvec.data[3];
    let gam1 = a * cth * k0 - a * a * cth * sth * sth * k3;
    let gam2 = r * (r * r + a * a) * sth * k3 - a * r * sth * k0;
    let gam3 = a * a * cth * sth * sth * k1 - r * (r * r + a * a) * sth * k2;
    let del1 = r * k0 - r * a * sth * sth * k3;
    // sign changed 8/24/2015 (upstream note; changes results by <1%)
    let del2 = -a * cth * sth * (r * r + a * a) * k3 + a * a * sth * cth * k0;
    let del3 = r * a * sth * sth * k1 + a * cth * sth * (r * r + a * a) * k2;
    let denom = (gam2 * del1 - gam1 * del2) * (g33 * k3 + g03 * k0)
        + (gam3 * del2 - gam2 * del3) * g11 * k1
        - (gam3 * del1 - gam1 * del3) * g22 * k2;
    let mut f1 = (gam2 * kap1 - del2 * kap2) * (g33 * k3 + g03 * k0)
        - g22 * k2 * (gam3 * kap1 - del3 * kap2);
    let mut f2 = (del1 * kap2 - gam1 * kap1) * (g33 * k3 + g03 * k0)
        + g11 * k1 * (gam3 * kap1 - del3 * kap2);
    let mut f3 = g22 * k2 * (gam1 * kap1 - del1 * kap2) - g11 * k1 * (gam2 * kap1 - del2 * kap2);
    if denom.abs() > 0.0 {
        f1 /= denom;
        f2 /= denom;
        f3 /= denom;
    }
    (f1, f2, f3)
}

/// Results of the comoving orthonormal frame construction.
#[derive(Clone, Copy, Debug, Default)]
pub struct ComovingOrtho {
    /// sin(2 xi): angle between the polarization basis and the magnetic
    /// field projection.
    pub s2xi: f64,
    /// cos(2 xi).
    pub c2xi: f64,
    /// angle between the wave vector and magnetic field.
    pub ang: f64,
    /// Doppler/gravitational redshift factor g = 1 / khat^t.
    pub g: f64,
    /// projection factor used by the thin-disk polarization model.
    pub cosne: f64,
}

/// Comoving orthonormal frame construction (upstream `comoving_ortho_core`;
/// Beckwith et al. 2008, Shcherbakov & Huang 2011).
///
/// `u`, `b`, `k` must carry the BL covariant metric as assigned by the fluid
/// model / ray construction (this function re-assigns them, matching
/// upstream's mutation).
#[allow(clippy::too_many_arguments)]
pub fn comoving_ortho(
    r: f64,
    th: f64,
    a: f64,
    alpha: f64,
    beta: f64,
    mus: f64,
    u: &mut FourVector,
    b: &mut FourVector,
    k: &mut FourVector,
) -> ComovingOrtho {
    let metric = blmetric_cov(r, th, a);
    let tmetric = metric; // transpose of a symmetric matrix
    let gtt = metric.0[0];
    let gtp = metric.0[3];
    let grr = metric.0[4];
    let gmm = metric.0[7];
    let gpp = metric.0[9];

    let ut = u.data[0];
    let ur = u.data[1];
    let um = u.data[2];
    let up = u.data[3];
    let utc = gtt * ut + gtp * up;
    let upc = gpp * up + gtp * ut;
    let urc = grr * ur;
    let umc = um * gmm;

    // Walker-Penrose transported basis
    let kap1 = alpha + a * (1.0 - mus * mus).sqrt();
    let kap2 = -beta;
    let (al1, al2, al3);
    if kap1 != 0.0 || kap2 != 0.0 {
        let (f1, f2, f3) = transport_perpk(k, r, th, a, &metric, kap1, kap2);
        al1 = f1;
        al2 = f2;
        al3 = f3;
    } else {
        al1 = 0.0;
        al2 = 0.0;
        al3 = 1.0 / metric.0[9].sqrt();
    }
    let mut aa = FourVector::new([0.0, al1, al2, al3], tmetric);

    // Kulkarni et al. (2011) normalizations
    let delta = r * r + a * a - 2.0 * r;
    let nr2 = -grr * (utc * ut + upc * up) * (1.0 + umc * um);
    let nm2 = gmm * (1.0 + umc * um);
    let np2 = -(utc * ut + upc * up) * delta * th.sin() * th.sin();
    let mut ekt = u.scale(-1.0);
    let ekr = FourVector::new(
        [
            urc * ut / nr2.sqrt(),
            -(utc * ut + upc * up) / nr2.sqrt(),
            0.0,
            urc * up / nr2.sqrt(),
        ],
        tmetric,
    );
    let ekm = FourVector::new(
        [
            umc * ut / nm2.sqrt(),
            umc * ur / nm2.sqrt(),
            (1.0 + umc * um) / nm2.sqrt(),
            umc * up / nm2.sqrt(),
        ],
        tmetric,
    );
    let ekp = FourVector::new([upc / np2.sqrt(), 0.0, 0.0, -utc / np2.sqrt()], tmetric);
    ekt.assign_metric(tmetric);

    let bhatt = ekt.dot(b);
    let bhatr = ekr.dot(b);
    let bhatm = ekm.dot(b);
    let bhatp = ekp.dot(b);
    let khatt = ekt.dot(k);
    let khatr = ekr.dot(k);
    let khatm = ekm.dot(k);
    let khatp = ekp.dot(k);
    let ahatt = ekt.dot(&aa);
    let ahatr = ekr.dot(&aa);
    let ahatm = ekm.dot(&aa);
    let ahatp = ekp.dot(&aa);

    let bhat = [bhatt, bhatr, bhatm, bhatp];
    let khat = [khatt, khatr, khatm, khatp];
    let ahat = [ahatt, ahatr, ahatm, ahatp];
    let _ = ahat;

    // upstream re-assigns the metric on the caller's vectors here
    b.assign_metric(tmetric);
    u.assign_metric(tmetric);
    let _ = &mut aa;

    let bdotb = b.dot(b);
    // upstream re-computes bdotk from the *frame* components (spatial dot
    // only; grtrans line 694), not as the BL inner product
    let bdotk = bhat[1] * khat[1] + bhat[2] * khat[2] + bhat[3] * khat[3];
    let om = -k.dot(u);
    let om2 = om * om;

    let aahat = [
        ahatr - khat[1] * ahatt / khat[0],
        ahatm - khat[2] * ahatt / khat[0],
        ahatp - khat[3] * ahatt / khat[0],
    ];
    let knorm = khatr * khatr + khatm * khatm + khatp * khatp;
    let sq = knorm.sqrt();
    let bbhat = [
        -(aahat[1] * khatp - aahat[2] * khatm) / sq,
        -(aahat[2] * khatr - aahat[0] * khatp) / sq,
        -(aahat[0] * khatm - aahat[1] * khatr) / sq,
    ];
    let _ = om2;

    let (s2xi, c2xi, angnorm);
    if bdotb > 0.0 {
        let aadotbp = bhat[1] * aahat[0] + bhat[2] * aahat[1] + bhat[3] * aahat[2];
        let bpdotbb = bhat[1] * bbhat[0] + bhat[2] * bbhat[1] + bhat[3] * bbhat[2];
        s2xi = -2.0 * aadotbp * bpdotbb / (aadotbp.powi(2) + bpdotbb.powi(2));
        c2xi = (bpdotbb * bpdotbb - aadotbp * aadotbp) / (aadotbp.powi(2) + bpdotbb.powi(2));
        angnorm = bdotk / knorm.sqrt() / bdotb.sqrt();
    } else {
        s2xi = 0.0;
        c2xi = 1.0;
        angnorm = 0.5;
    }
    let angnorm = angnorm.clamp(-0.99, 0.99);
    let ang = angnorm.acos();
    // g = 1 / khat^t in the comoving orthonormal frame
    let g = 1.0 / khat[0];
    let cosne = g * (beta * beta + mus * mus * (alpha * alpha - a * a)).sqrt() / r;

    ComovingOrtho {
        s2xi,
        c2xi,
        ang,
        g,
        cosne,
    }
}

/// Polarization four-vector for a thin disk in the LNRF (upstream
/// `calc_polvec`, from Agol 1997), with `mu = cos(theta)`.
pub fn calc_polvec(r: f64, mu: f64, p: &FourVector, a: f64, psi: f64) -> FourVector {
    let delta = r * r - 2.0 * r + a * a;
    let ar = (r * r + a * a).powi(2) - a * a * delta * (1.0 - mu * mu);
    let om = 2.0 * a * r / ar;
    let rho = r * r + a * a * mu * mu;
    let ptt = r * (delta / ar).sqrt() * p.data[0];
    let prt = r / delta.sqrt() * p.data[1];
    let ptht = r * p.data[2];
    let ppht = ar.sqrt() / r * (p.data[3] - om * p.data[0]);
    // LNRF velocity for a Keplerian disk
    let mut vel = 1.0 / (r.powf(1.5) + a);
    let epsi = (1.0 - mu * mu).sqrt() * (ar / rho).sqrt();
    let enu = (delta * rho / ar).sqrt();
    vel = epsi / enu * (vel - om);
    let frl = delta.sqrt() / r * (vel * (ptt - prt * prt / ptt) - ppht);
    let fthl = -vel * prt * ptht / ptt / r;
    let fphl = r * prt / ar.sqrt() * (1.0 - vel * ppht / ptt);
    let frp = delta.sqrt() * ptht * prt / r * (-1.0 + vel * ppht / ptt);
    let fthp = 1.0 / r
        * (prt * prt + (1.0 + vel * vel) * ppht * ppht - 2.0 * vel * ppht * ptt
            + vel * ptht * ptht * ppht / ptt);
    let fphp =
        r * ptht / ar.sqrt() * (-(1.0 + vel * vel) * ppht + vel * ptt + vel * ppht * ppht / ptt);
    let fr = psi.cos() * frl + psi.sin() * frp;
    let fth = psi.cos() * fthl + psi.sin() * fthp;
    let fph = psi.cos() * fphl + psi.sin() * fphp;
    let mut fourf = FourVector::flat([0.0, fr, fth, fph]);
    let th = mu.acos();
    fourf.assign_metric(blmetric_cov(r, th, a));
    let normf = fourf.norm_sq();
    for i in 0..4 {
        fourf.data[i] /= normf.sqrt();
    }
    fourf
}

/// Complex Penrose-Walker constant (upstream `calc_kappapw`), returned as
/// `(real, imag)`.
pub fn calc_kappapw(a: f64, r: f64, mu: f64, p: &FourVector, f: &FourVector) -> (f64, f64) {
    let alpha = (p.data[0] * f.data[1] - p.data[1] * f.data[0])
        + a * (1.0 - mu * mu) * (p.data[1] * f.data[3] - p.data[3] * f.data[1]);
    let beta =
        (r * r + a * a) * (1.0 - mu * mu).sqrt() * (p.data[3] * f.data[2] - p.data[2] * f.data[3])
            - a * (1.0 - mu * mu).sqrt() * (p.data[0] * f.data[2] - p.data[2] * f.data[0]);
    // kappapw = cmplx(alpha,-beta) * cmplx(r,-a*mu)
    let (re1, im1) = (alpha, -beta);
    let (re2, im2) = (r, -a * mu);
    (re1 * re2 - im1 * im2, re1 * im2 + im1 * re2)
}

/// Polarization angle and related quantities for single-point rays
/// (upstream `calc_polar_psi`). Returns `(s2psi, c2psi, cosne)`.
#[allow(clippy::too_many_arguments)]
pub fn calc_polar_psi(
    r: f64,
    muf: f64,
    q2: f64,
    a: f64,
    alpha: f64,
    beta: f64,
    rshift: f64,
    mus: f64,
    p: &FourVector,
) -> (f64, f64, f64) {
    let f = calc_polvec(r, muf, p, a, 0.0);
    let (kappa_re, kappa_im) = calc_kappapw(a, r, muf, p, &f);
    let kappa2 = kappa_re;
    let kappa1 = -kappa_im;
    let gammac = -alpha - a * (1.0 - mus * mus);
    let denom = beta * kappa2 - gammac * kappa1;
    let num = -beta * kappa1 - gammac * kappa2;
    let polarpsi = denom.atan2(num);
    let s2psi = (2.0 * polarpsi).sin();
    let c2psi = (2.0 * polarpsi).cos();
    let cosne = rshift * q2.sqrt() / r;
    (s2psi, c2psi, cosne)
}

/// `u^0` from three-velocities in BL coordinates (upstream `calc_u0`).
pub fn calc_u0(metric: &Metric, vr: f64, vth: f64, vph: f64) -> f64 {
    let m = &metric.0;
    (-1.0 / (m[0] + m[4] * vr * vr + m[7] * vth * vth + m[9] * vph * vph + 2.0 * m[3] * vph)).sqrt()
}

/// Four-velocity of a plunging geodesic with ISCO energy/angular momentum
/// (upstream `calc_plunging_vel`).
pub fn calc_plunging_vel(a: f64, r: f64) -> FourVector {
    let (ems, lms, _rms) = calc_rms_constants(a);
    let th = PI / 2.0;
    let metriccon = blmetric_con(r, th, a);
    let m = &metriccon.0;
    let pt = -m[0] * ems + m[3] * lms;
    let denom = -m[4] * (1.0 + m[0] * ems * ems - 2.0 * m[3] * ems * lms + m[9] * lms * lms);
    let pr = if denom > 0.0 { -denom.sqrt() } else { 0.0 };
    let pphi = -m[3] * ems + m[9] * lms;
    let mut fu = FourVector::flat([pt, pr, 0.0, pphi]);
    fu.assign_metric(metriccon);
    fu
}

/// Circular-orbit four-velocity at arbitrary theta, built from the
/// equatorial plunging solution (upstream `rms_vel`).
pub fn rms_vel(a: f64, th: f64, r: f64) -> FourVector {
    let fueq = calc_plunging_vel(a, r);
    let theq = PI / 2.0;
    let (vrl, vtl, vpl) = lnrf_frame(
        fueq.data[1] / fueq.data[0],
        fueq.data[2] / fueq.data[0],
        fueq.data[3] / fueq.data[0],
        r,
        a,
        theq,
    );
    let (vrr, vtt, vpp) = lnrf_frame_inv(vrl, vtl, vpl, r, a, th);
    let metric = blmetric_cov(r, th, a);
    let u0 = calc_u0(&metric, vrr, vtt, vpp);
    let mut fu = FourVector::new([u0, u0 * vrr, u0 * vtt, u0 * vpp], metric);
    fu.assign_metric(metric);
    fu
}

/// Radial surface integral of a 3-D field, assuming uniform phi
/// (upstream `surf_integral`).
///
/// Arrays have shape `(nx1, nx2, nx3)` in row-major order
/// (`x[i][j][k] = flat[(i*nx2 + j)*nx3 + k]`). Returns length `nx1`.
pub fn surf_integral(
    x: &[f64],
    r: &[f64],
    th: &[f64],
    ph: &[f64],
    a: f64,
    dims: (usize, usize, usize),
) -> Vec<f64> {
    let (nx1, nx2, nx3) = dims;
    assert_eq!(x.len(), nx1 * nx2 * nx3);
    let dph = ph[1] - ph[0];
    let idx = |i: usize, j: usize, k: usize| (i * nx2 + j) * nx3 + k;
    let mut s = vec![0.0f64; nx1];
    for i in 0..nx1 {
        let mut acc = 0.0f64;
        for j in 0..nx2 {
            // upstream: dth(:,2:nx2,:) = forward differences;
            // dth(:,1,:) = dth(:,nx2,:) (wraparound of the last difference)
            let dth = if j + 1 < nx2 {
                th[idx(i, j + 1, 0)] - th[idx(i, j, 0)]
            } else {
                th[idx(i, nx2 - 1, 0)] - th[idx(i, nx2 - 2, 0)]
            };
            for k in 0..nx3 {
                let rr = r[idx(i, j, k)];
                let tt = th[idx(i, j, k)];
                let gdet = (rr * rr + a * a * tt.cos() * tt.cos()) * tt.sin();
                acc += x[idx(i, j, k)] * gdet * dth * dph;
            }
        }
        s[i] = acc;
    }
    s
}

/// theta-integral variant (upstream `th_integral`): returns shape
/// `(nx1, nx3)` flattened row-major.
pub fn th_integral(
    x: &[f64],
    r: &[f64],
    th: &[f64],
    ph: &[f64],
    a: f64,
    dims: (usize, usize, usize),
) -> (Vec<f64>, usize, usize) {
    let (nx1, nx2, nx3) = dims;
    let dph = ph[1] - ph[0];
    let idx = |i: usize, j: usize, k: usize| (i * nx2 + j) * nx3 + k;
    let mut s = vec![0.0f64; nx1 * nx3];
    for i in 0..nx1 {
        for k in 0..nx3 {
            let mut acc = 0.0f64;
            for j in 0..nx2 {
                let dth = if j + 1 < nx2 {
                    th[idx(i, j + 1, 0)] - th[idx(i, j, 0)]
                } else {
                    th[idx(i, nx2 - 1, 0)] - th[idx(i, nx2 - 2, 0)]
                };
                let rr = r[idx(i, j, k)];
                let tt = th[idx(i, j, k)];
                let gdet = (rr * rr + a * a * tt.cos() * tt.cos()) * tt.sin();
                acc += x[idx(i, j, k)] * gdet * dth * dph;
            }
            s[i * nx3 + k] = acc * nx3 as f64;
        }
    }
    (s, nx1, nx3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_space_limit() {
        // far from the hole, BL -> Minkowski (corrections are O(1/r))
        let m = blmetric_cov(1e8, 1.0, 0.0);
        assert!((m.0[0] + 1.0).abs() < 1e-7);
        assert!((m.0[4] - 1.0).abs() < 1e-7);
    }

    #[test]
    fn metric_inverse() {
        // cov * con should be the identity for several points
        let a = 0.9;
        for &(r, th) in &[(3.0, 0.9), (10.0, 1.3), (100.0, 0.4)] {
            let cov = blmetric_cov(r, th, a).matrix();
            let con = blmetric_con(r, th, a).matrix();
            for i in 0..4 {
                for j in 0..4 {
                    let mut s = 0.0;
                    for k in 0..4 {
                        s += cov[i][k] * con[k][j];
                    }
                    let expect = if i == j { 1.0 } else { 0.0 };
                    assert!((s - expect).abs() < 1e-8, "mismatch at ({i},{j}): {s}");
                }
            }
        }
    }

    #[test]
    fn isco_schwarzschild() {
        // r_isco = 6M for a = 0
        let r = calc_rms(0.0);
        assert!((r - 6.0).abs() < 1e-3, "r_isco = {r}");
    }

    #[test]
    fn null_vector_is_null() {
        // a photon wave vector constructed by calc_nullp must be null
        let a = 0.9;
        let r = 10.0;
        let th: f64 = 1.2;
        let mu = th.cos();
        let (q2, l) = (5.0, 2.0);
        let kv = calc_nullp(q2, l, a, r, mu, 1.0, 1.0, true, true);
        let mut k = FourVector::new(kv, blmetric_cov(r, th, a));
        // raise and lower: k_mu k^mu with g^{mu nu} k_mu k_nu
        let con = blmetric_con(r, th, a);
        let mut kup = [0.0f64; 4];
        let m = &con.0;
        let kl = k.lower();
        kup[0] = m[0] * kl[0] + m[1] * kl[1] + m[2] * kl[2] + m[3] * kl[3];
        kup[1] = m[1] * kl[0] + m[4] * kl[1] + m[5] * kl[2] + m[6] * kl[3];
        kup[2] = m[2] * kl[0] + m[5] * kl[1] + m[7] * kl[2] + m[8] * kl[3];
        kup[3] = m[3] * kl[0] + m[6] * kl[1] + m[8] * kl[2] + m[9] * kl[3];
        let norm = kl[0] * kup[0] + kl[1] * kup[1] + kl[2] * kup[2] + kl[3] * kup[3];
        assert!(norm.abs() < 1e-9, "k^2 = {norm}");
        k.assign_metric(blmetric_cov(r, th, a));
    }

    #[test]
    fn comoving_frame_redshift_positive() {
        // A Keplerian disk flow should give a finite, positive g factor
        let a = 0.9;
        let r = 8.0;
        let th = PI_F64 / 2.0;
        let mut u = rms_vel(a, th, r);
        // magnetic field: radial-ish
        let mut b = FourVector::new([0.0, 1.0, 0.0, 0.0], blmetric_cov(r, th, a));
        // photon: use calc_nullp with arbitrary constants
        let kv = calc_nullp(4.0, 2.0, a, r, th.cos(), 1.0, 1.0, true, true);
        let mut k = FourVector::new(kv, blmetric_cov(r, th, a));
        let out = comoving_ortho(r, th, a, 3.0, 1.0, th.cos(), &mut u, &mut b, &mut k);
        assert!(out.g.is_finite() && out.g > 0.0, "g = {}", out.g);
        assert!(out.c2xi.abs() <= 1.0 + 1e-12);
        assert!(out.s2xi.abs() <= 1.0 + 1e-12);
    }
}
