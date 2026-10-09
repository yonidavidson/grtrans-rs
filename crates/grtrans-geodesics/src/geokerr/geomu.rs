//! GEOMU: given the final inverse radius, solve for the final mu and
//! classify the radial root case. Translation of `geokerr_wrapper.f`
//! GEOMU + FINDMROOTS + INDEP_MUF.
//!
//! Upstream's in/out argument convention is modelled with
//! [`GeomuInputs`] (cached first-point data: quartic roots, NCASE, H1 and
//! the cached integrals) and [`GeomuResult`].

use crate::geokerr::elliptic::{self, asech};
use crate::geokerr::imu;
use crate::geokerr::special::sncndn;
use num_complex::Complex64;
use std::f64::consts::PI;

/// Cached first-point data, passed back into GEOMU on subsequent calls
/// (upstream reads these arguments as inputs when FIRSTPT is false).
#[derive(Clone, Copy, Debug, Default)]
pub struct GeomuInputs {
    pub h1: f64,
    pub u1: f64,
    pub u2: f64,
    pub u3: f64,
    pub u4: f64,
    pub rffu0: f64,
    pub rffu1: f64,
    pub rffmu1: f64,
    pub rffmu2: f64,
    pub rffmu3: f64,
    pub iu0: f64,
    pub i1mu: f64,
    pub i3mu: f64,
    pub ncase: i32,
}

impl GeomuInputs {
    /// Initial cache: matches the uninitialized-input call pattern of the
    /// first GEOMU call (FIRSTPT=true), which overwrites everything used.
    pub fn first() -> Self {
        Self::default()
    }
}

/// Output of one GEOMU call (superset of upstream's out/in arguments).
#[derive(Clone, Copy, Debug, Default)]
pub struct GeomuResult {
    /// new value of the final mu (unchanged from input in some branches,
    /// exactly as upstream)
    pub muf: f64,
    pub iu: f64,
    pub tpm: i32,
    pub tpr: i32,
    pub h1: f64,
    pub u1: f64,
    pub u2: f64,
    pub u3: f64,
    pub u4: f64,
    pub rffu0: f64,
    pub rffu1: f64,
    pub rffmu1: f64,
    pub rffmu2: f64,
    pub rffmu3: f64,
    pub iu0: f64,
    pub i1mu: f64,
    pub i3mu: f64,
    pub ncase: i32,
    /// possibly-modified uf (upstream modifies UF in place when clamping)
    pub uf: f64,
    /// possibly-modified u0 (upstream modifies U0 in unphysical branches)
    pub u0: f64,
    /// first-point cache to pass into subsequent calls
    pub cache: GeomuInputs,
}

/// Fortran `(-1)**n` for integer n (parity; negative n handled too).
#[inline]
fn neg_one_pow(n: i32) -> f64 {
    if n % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

/// Fortran `INT(x)` as i32 (truncation toward zero).
#[inline]
fn f_int(x: f64) -> i32 {
    x as i32
}

/// Roots of M(mu) = 0 (upstream `findmroots`). Returns `(muminus, muplus)`.
pub fn findmroots(q2: f64, l: f64, a: f64, mu0: f64) -> (f64, f64) {
    let ql2 = q2 + l * l;
    let s = if a * a - ql2 >= 0.0 { 1.0 } else { -1.0 };
    let yy = -0.5 * (a * a - ql2 + s * ((a * a - ql2).powi(2) + 4.0 * q2 * a * a).sqrt());
    let (mneg, mut mpos) = if (a * a - ql2) < 0.0 {
        (-yy / a / a, q2 / yy)
    } else {
        (q2 / yy, -yy / a / a)
    };
    mpos = mpos.min(1.0);
    if mneg > 0.0 {
        // asymmetric roots
        let sgn = if mu0 >= 0.0 { 1.0 } else { -1.0 };
        let mut muplus = sgn * mpos.sqrt();
        let mut muminus = sgn * mneg.sqrt();
        if muminus.abs() > mu0.abs() {
            muminus = mu0;
        }
        if mu0.abs() > muplus.abs() {
            muplus = mu0;
        }
        (muminus, muplus)
    } else {
        // symmetric roots: orbit can cross the equatorial plane
        let mut muplus = mpos.sqrt();
        if muplus < mu0 {
            muplus = mu0;
        }
        (-muplus, muplus)
    }
}

/// Upstream `indep_muf`: mu value at index K between MU0 and MUF including
/// turning points. Returns `(mun, tpmk)`.
#[allow(clippy::too_many_arguments)]
pub fn indep_muf(
    muminus: f64,
    muplus: f64,
    mu0: f64,
    muf: f64,
    tpm0: i32,
    tpmf: i32,
    sm: f64,
    k: i32,
    npts: i32,
    offset: f64,
) -> (f64, i32) {
    let mun;
    let mut tpmk;
    if muf <= muplus && muf >= muminus {
        let dtpm = tpmf - tpm0;
        let a1 = sm * neg_one_pow(tpm0);
        let a2 = sm * neg_one_pow(tpmf);
        let a3 = 2.0 * ((2.0 * dtpm as f64 + 3.0 - a1) / 4.0).floor() - 1.0;
        let dmu = a1 * (muplus - mu0) + a2 * (muf - muminus) + a3 * (muplus - muminus);
        let mu1 = (a1 + 1.0) / 2.0 * muplus + (1.0 - a1) / 2.0 * muminus;
        let mu2 = (1.0 - a1) / 2.0 * muplus + (a1 + 1.0) / 2.0 * muminus;
        let dmu1 = (mu1 - mu0).abs();
        let dmu2 = (mu2 - mu1).abs();
        let dmuk = (k as f64 - offset) * dmu / npts as f64;
        let mut dmukdmu1 = f_int(dmuk / dmu1);
        if dmukdmu1 < 0 {
            dmukdmu1 = 1;
        }
        tpmk = dmukdmu1.min(1) + f_int((dmuk - dmu1) / dmu2).max(0);
        let a1k = sm * neg_one_pow(tpm0);
        let a2k = sm * neg_one_pow(tpmk + tpm0);
        let a3k = 2.0 * f_int((2.0 * tpmk as f64 + 3.0 - a1k) / 4.0) as f64 - 1.0;
        mun = muminus + 1.0 / a2k * (dmuk - a1k * (muplus - mu0) - a3k * (muplus - muminus));
    } else {
        // no solution exists
        mun = 0.0;
        tpmk = 0;
    }
    tpmk += tpm0;
    (mun, tpmk)
}

/// Upstream `geomu`. `cached` carries first-point outputs; on
/// `firstpt == true` they are recomputed and returned in the result.
#[allow(clippy::too_many_arguments)]
pub fn geomu(
    u0: f64,
    uf_in: f64,
    mu0: f64,
    muf_in: f64,
    a: f64,
    l: f64,
    l2: f64,
    q2: f64,
    tpm_in: i32,
    tpr_in: i32,
    su: f64,
    sm: f64,
    cached: &GeomuInputs,
    pht: bool,
    firstpt: bool,
) -> GeomuResult {
    let one = 1.0f64;
    let two = 2.0f64;
    let third = 0.3333333333333333f64;
    let pi = PI;
    let _pi2 = two * pi;
    let uplus = one / (one + (one - a * a).sqrt());
    let ql2 = q2 + l2;

    // Eq. (13): coefficients of the quartic in U(u)
    let cc = a * a - q2 - l2;
    let dd = two * ((a - l) * (a - l) + q2);
    let ee = -a * a * q2;

    let mut out = GeomuResult {
        muf: muf_in,
        uf: uf_in,
        u0,
        tpr: tpr_in,
        h1: cached.h1,
        u1: cached.u1,
        u2: cached.u2,
        u3: cached.u3,
        u4: cached.u4,
        rffu0: cached.rffu0,
        rffu1: cached.rffu1,
        rffmu1: cached.rffmu1,
        rffmu2: cached.rffmu2,
        rffmu3: cached.rffmu3,
        iu0: cached.iu0,
        i1mu: cached.i1mu,
        i3mu: cached.i3mu,
        ncase: cached.ncase,
        ..Default::default()
    };
    let mut u0 = u0;
    let mut uf = uf_in;
    let mut muf = muf_in;
    let iu;
    let mut iu1 = 0.0f64;
    #[allow(unused_assignments)]
    let mut muminus = 0.0f64;
    let mut i1mu = cached.i1mu;
    let mut i3mu = cached.i3mu;
    let mut rffmu1 = cached.rffmu1;
    let mut rffmu3 = cached.rffmu3;
    let mut rffmu2 = cached.rffmu2;
    let mut h1 = cached.h1;
    let mut tpm = tpm_in;

    if ee == 0.0 && dd != 0.0 {
        // ------------------------------------------------------------
        // cubic cases
        // ------------------------------------------------------------
        let p: [i32; 5] = [-1, -1, -1, 0, 0];
        let qq = cc * cc / dd / dd / 9.0;
        let rr = (two * cc.powi(3) / dd.powi(3) + 27.0 / dd) / 54.0;
        let dis = rr * rr - qq.powi(3);
        if dis < -1e-16 {
            // cubic real roots, u1 < 0 < u2 <= u3
            let theta = (rr / qq.powf(1.5)).acos();
            out.u1 = -two * qq.sqrt() * (theta / 3.0).cos() - cc / dd / 3.0;
            out.u2 = -two * qq.sqrt() * ((theta - two * pi) / 3.0).cos() - cc / dd / 3.0;
            out.u3 = -two * qq.sqrt() * ((theta + two * pi) / 3.0).cos() - cc / dd / 3.0;
            let mut u0_work = u0;
            let mut uf_work = uf;
            loop {
                // label 80 (retried when uf lands in the forbidden region)
                u0_work = u0;
                uf_work = uf;
                if u0_work <= out.u2 {
                    if uf_work > out.u2 {
                        uf_work = out.u2;
                    }
                    out.ncase = 1;
                    // Table 1 Row 1
                    if firstpt && u0_work != out.u2 {
                        let mut rff = out.rffu0;
                        out.iu0 = elliptic::ellcubicreal(
                            &p, -out.u1, one, out.u2, -one, out.u3, -one, 0.0, 0.0, &mut rff,
                            u0_work, out.u2,
                        );
                        out.rffu0 = rff;
                    } else if u0_work == out.u2 {
                        out.iu0 = 0.0;
                    }
                    if uf_work != out.u2 {
                        let mut rff = out.rffu1;
                        iu1 = elliptic::ellcubicreal(
                            &p, -out.u1, one, out.u2, -one, out.u3, -one, 0.0, 0.0, &mut rff,
                            uf_work, out.u2,
                        );
                        out.rffu1 = rff;
                    }
                    iu = su * (out.iu0 - neg_one_pow(out.tpr) * iu1) / dd.sqrt();
                    break;
                } else if u0_work >= out.u3 {
                    if uf_work < out.u3 {
                        uf_work = out.u3;
                    }
                    out.ncase = 2;
                    // Table 1 Row 2
                    if firstpt && u0_work != out.u3 {
                        let mut rff = out.rffu0;
                        out.iu0 = -elliptic::ellcubicreal(
                            &p, -out.u1, one, -out.u2, one, -out.u3, one, 0.0, 0.0, &mut rff,
                            out.u3, u0_work,
                        );
                        out.rffu0 = rff;
                    } else if u0_work == out.u3 {
                        out.iu0 = 0.0;
                    }
                    if uf_work != out.u3 {
                        let mut rff = out.rffu1;
                        iu1 = -elliptic::ellcubicreal(
                            &p, -out.u1, one, -out.u2, one, -out.u3, one, 0.0, 0.0, &mut rff,
                            out.u3, uf_work,
                        );
                        out.rffu1 = rff;
                    }
                    iu = su * (out.iu0 - neg_one_pow(out.tpr) * iu1) / dd.sqrt();
                    break;
                } else {
                    // uf in the forbidden region u2 < uf < u3: upstream
                    // issues a warning and modifies u0, then retries.
                    if su == 1.0 {
                        u0 = out.u3;
                    } else {
                        u0 = out.u2;
                    }
                    // loop retries (upstream GOTO 80)
                }
            }
            uf = uf_work;
        } else if dis.abs() < 1e-16 {
            // cubic with equal roots: elementary functions
            out.u1 = -two * qq.sqrt() - cc / dd / 3.0;
            out.u2 = -two * qq.sqrt() * (two * pi / 3.0).cos() - cc / dd / 3.0;
            out.u3 = out.u2;
            out.tpr = 0;
            if u0 <= out.u2 {
                out.ncase = 1;
                if uf > out.u2 {
                    uf = out.u2;
                    iu = su * 1e300;
                } else {
                    let farg = ((uf - out.u1).sqrt() + (out.u2 - out.u1).sqrt())
                        / ((uf - out.u1).sqrt() - (out.u2 - out.u1).sqrt()).abs();
                    let sarg = ((u0 - out.u1).sqrt() + (out.u2 - out.u1).sqrt())
                        / ((u0 - out.u1).sqrt() - (out.u2 - out.u1).sqrt()).abs();
                    iu = su * (farg.ln() - sarg.ln()) / ((out.u2 - out.u1) * dd).sqrt();
                }
            } else if u0 >= out.u2 {
                out.ncase = 2;
                if uf < out.u2 {
                    uf = out.u2;
                    iu = su * 1e300;
                } else {
                    let farg = ((uf - out.u1).sqrt() + (out.u2 - out.u1).sqrt())
                        / ((uf - out.u1).sqrt() - (out.u2 - out.u1).sqrt()).abs();
                    let sarg = ((u0 - out.u1).sqrt() + (out.u2 - out.u1).sqrt())
                        / ((u0 - out.u1).sqrt() - (out.u2 - out.u1).sqrt()).abs();
                    iu = -su * (farg.ln() - sarg.ln()) / ((out.u2 - out.u1) * dd).sqrt();
                }
            } else {
                iu = 0.0;
            }
        } else {
            // cubic complex case with one real root
            out.ncase = 3;
            let sign_rr = if rr >= 0.0 { 1.0 } else { -1.0 };
            let aa = -sign_rr * (rr.abs() + dis.sqrt()).powf(third);
            let bb = if aa != 0.0 { qq / aa } else { 0.0 };
            out.u1 = (aa + bb) - cc / dd / 3.0;
            let f = -one / dd / out.u1;
            let g = f / out.u1;
            if uf > u0 {
                let mut rff = out.rffu0;
                iu =
                    su * elliptic::ellcubiccomplex(
                        &p, -out.u1, one, 0.0, 0.0, f, g, one, &mut rff, u0, uf,
                    ) / dd.sqrt();
                out.rffu0 = rff;
            } else if uf < u0 {
                let mut rff = out.rffu0;
                iu =
                    -su * elliptic::ellcubiccomplex(
                        &p, -out.u1, one, 0.0, 0.0, f, g, one, &mut rff, uf, u0,
                    ) / dd.sqrt();
                out.rffu0 = rff;
            } else {
                iu = 0.0;
            }
        }
        // muf for the cubic cases (only possible when a=0 or q2=0)
        if q2 == 0.0 {
            let s1 = if mu0 >= 0.0 { 1.0 } else { -1.0 };
            if l2 < a * a && mu0 != 0.0 {
                let muplus = s1 * (one - l2 / a / a).sqrt();
                muf = muplus / ((a * muplus).abs() * iu - s1 * sm * asech(mu0 / muplus)).cosh();
            } else {
                muf = 0.0;
            }
        } else if a == 0.0 {
            let a1 = sm;
            let mut muplus = (q2 / ql2).sqrt();
            if mu0 > muplus {
                muplus = mu0;
            }
            i1mu = (mu0 / muplus).acos() / ql2.sqrt();
            i3mu = pi / ql2.sqrt();
            if sm == 1.0 {
                let ratio = (iu - i1mu) / i3mu;
                tpm = f_int(ratio) + f_int((1.0 + ratio.signum()) / 2.0);
            } else {
                tpm = f_int((iu + i1mu) / i3mu);
            }
            let a2 = sm * neg_one_pow(tpm);
            let a3 = two * f_int((two * tpm as f64 + 3.0 - sm) / 4.0) as f64 - one;
            muf = -muplus * (ql2.sqrt() * (iu - a1 * i1mu - a3 * i3mu) / a2).cos();
        }
    } else if ee == 0.0 && dd == 0.0 {
        // special case q2 = 0 and l = a: U(u) = 1
        out.ncase = 4;
        out.tpr = 0;
        iu = su * (uf - u0);
        muf = 0.0;
    } else {
        // ------------------------------------------------------------
        // quartic cases
        // ------------------------------------------------------------
        let p: [i32; 5] = [-1, -1, -1, -1, 0];
        if firstpt {
            let c: [Complex64; 5] = [
                Complex64::new(one, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(cc, 0.0),
                Complex64::new(dd, 0.0),
                Complex64::new(ee, 0.0),
            ];
            let root = crate::geokerr::special::zroots(&c, 4, true);
            let mut nreal = 0;
            for r in root.iter() {
                if r.im == 0.0 {
                    nreal += 1;
                }
            }
            if nreal == 2 {
                out.ncase = 5;
                out.u1 = root[0].re;
                if root[1].im == 0.0 {
                    out.u4 = root[1].re;
                } else {
                    out.u4 = root[3].re;
                }
            } else if nreal == 0 {
                out.ncase = 6;
            } else {
                out.u1 = root[0].re;
                out.u2 = root[1].re;
                out.u3 = root[2].re;
                out.u4 = root[3].re;
                loop {
                    // label 90
                    if out.u2 > uplus && out.u3 > uplus {
                        out.ncase = 5;
                        break;
                    } else if u0 <= out.u2 {
                        out.ncase = 7;
                        break;
                    } else if u0 >= out.u3 {
                        out.ncase = 8;
                        break;
                    } else {
                        // unphysical quartic real: modify u0 and retry
                        if su == 1.0 {
                            u0 = 0.0;
                        } else {
                            u0 = uplus;
                        }
                    }
                }
            }
        }
        if out.ncase == 5 {
            // quartic complex with 2 real roots; Table 1 Row 5
            let qs = if q2 >= 0.0 { 1.0 } else { -1.0 };
            let f = -qs * one / ee.abs() / out.u1 / out.u4;
            let g = (out.u4 + out.u1) / out.u1 / out.u4 * f;
            if uf > u0 {
                let mut rff = out.rffu0;
                iu =
                    su * elliptic::ellquarticcomplex(
                        &p,
                        -out.u1,
                        one,
                        qs * out.u4,
                        -qs * one,
                        0.0,
                        0.0,
                        f,
                        g,
                        one,
                        &mut rff,
                        u0,
                        uf,
                    ) / ee.abs().sqrt();
                out.rffu0 = rff;
            } else if uf < u0 {
                let mut rff = out.rffu0;
                iu =
                    -su * elliptic::ellquarticcomplex(
                        &p,
                        -out.u1,
                        one,
                        qs * out.u4,
                        -qs * one,
                        0.0,
                        0.0,
                        f,
                        g,
                        one,
                        &mut rff,
                        uf,
                        u0,
                    ) / ee.abs().sqrt();
                out.rffu0 = rff;
            } else {
                iu = 0.0;
            }
        } else if out.ncase == 6 {
            // quartic complex with no real roots; Table 1 Row 6
            let coefs: [Complex64; 7] = [
                Complex64::new(one, 0.0),
                Complex64::new(-cc / ee.sqrt(), 0.0),
                Complex64::new(-one, 0.0),
                Complex64::new(ee.sqrt() * (two * cc / ee - (dd / ee).powi(2)), 0.0),
                Complex64::new(-one, 0.0),
                Complex64::new(-cc / ee.sqrt(), 0.0),
                Complex64::new(one, 0.0),
            ];
            let hroots = crate::geokerr::special::zroots(&coefs, 6, true);
            h1 = 0.0;
            for r in hroots.iter() {
                if r.im == 0.0 {
                    h1 = r.re;
                    break;
                }
            }
            // upstream loops until a real root is found and then stops on
            // the first nonzero value; replicate the scan
            if h1 == 0.0 {
                for r in hroots.iter() {
                    if r.im == 0.0 && r.re != 0.0 {
                        h1 = r.re;
                        break;
                    }
                }
            }
            let h2 = one / h1;
            let g1 = dd / ee / (h2 - h1);
            let g2 = -g1;
            let f1 = one / ee.sqrt();
            let f2 = f1;
            if uf > u0 {
                let mut rff = out.rffu0;
                iu =
                    su * elliptic::elldoublecomplex(
                        &p, f1, g1, h1, f2, g2, h2, 0.0, 0.0, &mut rff, u0, uf,
                    ) / ee.abs().sqrt();
                out.rffu0 = rff;
            } else if uf < u0 {
                let mut rff = out.rffu0;
                iu =
                    -su * elliptic::elldoublecomplex(
                        &p, f1, g1, h1, f2, g2, h2, 0.0, 0.0, &mut rff, uf, u0,
                    ) / ee.abs().sqrt();
                out.rffu0 = rff;
            } else {
                iu = 0.0;
            }
        } else {
            // quartic real roots
            if (out.u3 - out.u2).abs() > 1e-12 {
                iu1 = 0.0;
                if out.ncase == 7 {
                    // Table 1 Row 7
                    if uf > out.u2 {
                        uf = out.u2;
                    }
                    if firstpt && u0 != out.u2 {
                        let mut rff = out.rffu0;
                        let iu0v = elliptic::ellquarticreal(
                            &p, -out.u1, one, out.u2, -one, out.u3, -one, out.u4, -one, 0.0, 0.0,
                            &mut rff, u0, out.u2,
                        );
                        out.iu0 = iu0v;
                        out.rffu0 = rff;
                    } else if u0 == out.u2 {
                        out.iu0 = 0.0;
                    }
                    if uf != out.u2 {
                        let mut rff = out.rffu1;
                        iu1 = elliptic::ellquarticreal(
                            &p, -out.u1, one, out.u2, -one, out.u3, -one, out.u4, -one, 0.0, 0.0,
                            &mut rff, uf, out.u2,
                        );
                        out.rffu1 = rff;
                    }
                    iu = su * (out.iu0 - neg_one_pow(out.tpr) * iu1) / ee.abs().sqrt();
                } else if out.ncase == 8 {
                    // Table 1 Row 8
                    if uf < out.u3 {
                        uf = out.u3;
                    }
                    if firstpt && u0 != out.u3 {
                        let mut rff = out.rffu0;
                        let iu0v = -elliptic::ellquarticreal(
                            &p, -out.u1, one, -out.u2, one, -out.u3, one, out.u4, -one, 0.0, 0.0,
                            &mut rff, out.u3, u0,
                        );
                        out.iu0 = iu0v;
                        out.rffu0 = rff;
                    } else if u0 == out.u3 {
                        out.iu0 = 0.0;
                    }
                    if uf != out.u3 {
                        let mut rff = out.rffu1;
                        iu1 = -elliptic::ellquarticreal(
                            &p, -out.u1, one, -out.u2, one, -out.u3, one, out.u4, -one, 0.0, 0.0,
                            &mut rff, out.u3, uf,
                        );
                        out.rffu1 = rff;
                    }
                    iu = su * (out.iu0 - neg_one_pow(out.tpr) * iu1) / ee.abs().sqrt();
                } else {
                    iu = 0.0;
                }
            } else {
                // equal roots quartic cases
                out.tpr = 0;
                if out.ncase == 7 {
                    if uf > out.u2 {
                        uf = out.u2;
                    }
                    let mut rff = out.rffu0;
                    out.iu0 = elliptic::ellquarticreal(
                        &p, -out.u1, one, out.u2, -one, out.u3, -one, out.u4, -one, 0.0, 0.0,
                        &mut rff, u0, uf,
                    );
                    out.rffu0 = rff;
                    iu = su * out.iu0 / ee.abs().sqrt();
                } else if out.ncase == 8 {
                    if uf < out.u3 {
                        uf = out.u3;
                    }
                    let mut rff = out.rffu0;
                    out.iu0 = -elliptic::ellquarticreal(
                        &p, -out.u1, one, -out.u2, one, -out.u3, one, out.u4, -one, 0.0, 0.0,
                        &mut rff, u0, uf,
                    );
                    out.rffu0 = rff;
                    iu = su * out.iu0 / ee.abs().sqrt();
                } else {
                    out.ncase = 0;
                    iu = 0.0;
                }
            }
        }
        // roots of M(mu) (Eq. (12))
        let yy = -0.5
            * (a * a - ql2
                + if a * a - ql2 >= 0.0 { 1.0 } else { -1.0 }
                    * ((a * a - ql2).powi(2) + 4.0 * q2 * a * a).sqrt());
        let (mut mneg, mut mpos) = if (a * a - ql2) < 0.0 {
            (-yy / a / a, q2 / yy)
        } else {
            (q2 / yy, -yy / a / a)
        };
        if mpos > 1.0 {
            mpos = 1.0;
        }
        let mut muplus = mpos.sqrt();
        if mneg < 0.0 {
            // symmetric roots case
            if muplus < mu0 {
                muplus = mu0;
                mpos = muplus * muplus;
            }
            muminus = -muplus;
            if firstpt {
                let vals = imu::calcimusym(a, mneg, mpos, mu0, muplus);
                i1mu = vals.imu0;
                i3mu = vals.imum;
                rffmu1 = vals.rf0;
                rffmu3 = vals.rfc;
            }
            let a1 = sm;
            if sm == 1.0 {
                let ratio = (iu - i1mu) / i3mu;
                tpm = f_int(ratio) + f_int((1.0 + ratio.signum()) / 2.0);
            } else {
                tpm = f_int((iu + i1mu) / i3mu);
            }
            let a2 = sm * neg_one_pow(tpm);
            let a3 = two * f_int((two * tpm as f64 + 3.0 - sm) / 4.0) as f64 - one;
            let iarg = a.abs() / a2 * (iu - a1 * i1mu - a3 * i3mu);
            let uarg = (mpos - mneg).sqrt() * iarg;
            let m1 = -mneg / (mpos - mneg);
            let (_sn, cn, _dn) = sncndn(uarg, m1);
            muf = muminus * cn;
            if pht {
                let vals = imu::calcimusymf(a, mneg, mpos, muf, muplus, i3mu);
                let _ = vals.imu0;
                i3mu = vals.imum;
                rffmu2 = vals.rf0;
            }
        } else {
            // asymmetric roots case
            let s1 = if mu0 >= 0.0 { 1.0 } else { -1.0 };
            muplus *= s1;
            muminus = s1 * mneg.sqrt();
            if muminus.abs() > mu0.abs() {
                muminus = mu0;
                mneg = muminus * muminus;
            }
            if mu0.abs() > muplus.abs() {
                muplus = mu0;
                mpos = muplus * muplus;
            }
            if firstpt {
                let vals = imu::calcimuasym(a, mneg, mpos, mu0, muplus);
                i1mu = vals.imu0;
                i3mu = vals.imum;
                rffmu1 = vals.rf0;
                rffmu3 = vals.rfc;
            }
            let a1 = sm;
            if sm == 1.0 {
                let ratio = (iu - i1mu) / i3mu;
                tpm = f_int(ratio) + f_int(((1.0 + ratio.signum()) / 2.0).abs());
            } else {
                tpm = f_int((iu + i1mu) / i3mu);
            }
            let a2 = sm * neg_one_pow(tpm);
            let a3 = two * f_int((two * tpm as f64 + 3.0 - sm) / 4.0) as f64 - one;
            let iarg = a.abs() / a2 * (iu - a1 * i1mu - a3 * i3mu);
            let uarg = muplus.abs() * iarg;
            let m1 = mneg / mpos;
            let (_sn, _cn, dn) = sncndn(uarg, m1);
            muf = muminus / dn;
            if pht {
                let vals = imu::calcimuasymf(a, mneg, mpos, muf, muplus, i3mu);
                let _ = vals.imu0;
                i3mu = vals.imum;
                rffmu2 = vals.rf0;
            }
        }
    }

    out.iu = iu;
    out.muf = muf;
    out.tpm = tpm;
    out.uf = uf;
    out.u0 = u0;
    out.i1mu = i1mu;
    out.i3mu = i3mu;
    out.rffmu1 = rffmu1;
    out.rffmu2 = rffmu2;
    out.rffmu3 = rffmu3;
    out.h1 = h1;
    out.cache = GeomuInputs {
        h1: out.h1,
        u1: out.u1,
        u2: out.u2,
        u3: out.u3,
        u4: out.u4,
        rffu0: out.rffu0,
        rffu1: out.rffu1,
        rffmu1: out.rffmu1,
        rffmu2: out.rffmu2,
        rffmu3: out.rffmu3,
        iu0: out.iu0,
        i1mu: out.i1mu,
        i3mu: out.i3mu,
        ncase: out.ncase,
    };
    out
}
