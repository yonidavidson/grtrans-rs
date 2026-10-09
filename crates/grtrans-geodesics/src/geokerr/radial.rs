//! GEOR: solve for the final inverse radius UF given the final polar angle
//! MUF, including root classification and the radial integrals.
//!
//! Direct translation of `geokerr_wrapper.f` subroutine GEOR
//! (lines 2672-3134).
//!
//! Naming note: upstream GEOR's output argument is named `IMU` (the
//! mu-integral used by GEOPHITIME as its `IU` argument); the local radial
//! integral `IU` computed inside GEOR is not returned. The Rust field
//! [`GeorResult::iu`] carries the IMU output so that the caller's variable
//! naming matches GEOKERR.

use crate::geokerr::elliptic::{self, asech};
use crate::geokerr::imu;
use crate::geokerr::special::sncndn;
use num_complex::Complex64;
use std::f64::consts::PI;

/// Cached first-point data passed back into GEOR on subsequent calls.
#[derive(Clone, Copy, Debug, Default)]
pub struct GeorInputs {
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
    pub i2mu: f64,
    pub i3mu: f64,
    pub ncase: i32,
}

/// Result of a GEOR call.
#[derive(Clone, Copy, Debug, Default)]
pub struct GeorResult {
    /// GEOR's IMU output (GEOPHITIME's IU)
    pub iu: f64,
    /// possibly modified final inverse radius (negative marks an invalid ray)
    pub uf: f64,
    /// possibly modified number of u turning points
    pub tpr: i32,
    pub ncase: i32,
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
    pub i2mu: f64,
    pub i3mu: f64,
}

impl GeorResult {
    /// Cached first-point data for subsequent GEOR calls.
    pub fn cache(&self) -> GeorInputs {
        GeorInputs {
            h1: self.h1,
            u1: self.u1,
            u2: self.u2,
            u3: self.u3,
            u4: self.u4,
            rffu0: self.rffu0,
            rffu1: self.rffu1,
            rffmu1: self.rffmu1,
            rffmu2: self.rffmu2,
            rffmu3: self.rffmu3,
            iu0: self.iu0,
            i1mu: self.i1mu,
            i2mu: self.i2mu,
            i3mu: self.i3mu,
            ncase: self.ncase,
        }
    }
}

#[inline]
fn fsign(x: f64) -> f64 {
    if x >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

#[inline]
fn f_int(x: f64) -> i32 {
    x as i32
}

/// Fortran `(-1)**n` for integer n.
#[inline]
fn neg_one_pow(n: i32) -> f64 {
    if n % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

/// Upstream `geor`.
#[allow(clippy::too_many_arguments)]
pub fn geor(
    u0: f64,
    uf_in: f64,
    mu0: f64,
    muf: f64,
    a: f64,
    l: f64,
    l2: f64,
    q2: f64,
    tpm: i32,
    tpr_in: i32,
    su: f64,
    sm: f64,
    cached: &GeorInputs,
    pht: bool,
    firstpt: bool,
) -> GeorResult {
    let one = 1.0f64;
    let two = 2.0f64;
    let half = 0.5f64;
    let third = 0.3333333333333333f64;
    let three = 3.0f64;
    let pi = PI;
    let uplus = one / (one + (one - a * a).sqrt());
    let ql2 = q2 + l2;
    let cc = a * a - q2 - l2;
    let dd = two * ((a - l) * (a - l) + q2);
    let ee = -a * a * q2;

    let mut out = GeorResult {
        uf: uf_in,
        tpr: tpr_in,
        ncase: cached.ncase,
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
        i2mu: cached.i2mu,
        i3mu: cached.i3mu,
        ..Default::default()
    };
    let mut u0 = u0;
    let mut uf = uf_in;
    let mut tpr = tpr_in;
    let mut imu = 0.0f64;
    let mut i1mu = cached.i1mu;
    let mut i2mu = cached.i2mu;
    let mut i3mu = cached.i3mu;
    let _iu0 = cached.iu0;

    if q2 == 0.0 {
        // Calculate imu in the special case q2=0: 0 or 1 mu turning points.
        if l2 >= a * a || mu0 == 0.0 {
            out.uf = -1.0;
            out.ncase = 0;
            out.tpr = tpr;
            out.iu = 0.0;
            return out;
        } else {
            let s1 = fsign(mu0);
            let a1 = s1 * sm;
            let a2 = s1 * sm * neg_one_pow(tpm + 1);
            let muplus = s1 * (one - l2 / a / a).sqrt();
            i1mu = one / (a * muplus).abs() * asech(mu0 / muplus);
            i2mu = one / (a * muplus).abs() * asech(muf / muplus);
            imu = a1 * i1mu + a2 * i2mu;
        }
    } else if a == 0.0 {
        // Calculate imu in the special case a=0
        let a1 = sm;
        let mut muplus = (q2 / ql2).sqrt();
        if mu0 > muplus {
            muplus = mu0;
        }
        i1mu = (half * pi - (mu0 / muplus).asin()) / ql2.sqrt();
        i3mu = pi / ql2.sqrt();
        let a2 = sm * neg_one_pow(tpm);
        let a3 = two * f_int((two * tpm as f64 + 3.0 - sm) / 4.0) as f64 - one;
        i2mu = (half * pi + (muf / muplus).asin()) / ql2.sqrt();
        imu = a1 * i1mu + a2 * i2mu + a3 * i3mu;
    }

    if ee == 0.0 && dd != 0.0 {
        // ----------------------------------------------------------------
        // cubic cases
        // ----------------------------------------------------------------
        let p: [i32; 5] = [-1, -1, -1, 0, 0];
        let qq = cc * cc / dd / dd / 9.0;
        let rr = (two * cc.powi(3) / dd.powi(3) + 27.0 / dd) / 54.0;
        let dis = rr * rr - qq.powi(3);
        if dis < -1e-16 {
            // cubic real roots with u1 < 0 < u2 <= u3
            let theta = (rr / qq.powf(1.5)).acos();
            out.u1 = -two * qq.sqrt() * (theta / 3.0).cos() - cc / dd / 3.0;
            out.u2 = -two * qq.sqrt() * ((theta - two * pi) / 3.0).cos() - cc / dd / 3.0;
            out.u3 = -two * qq.sqrt() * ((theta + two * pi) / 3.0).cos() - cc / dd / 3.0;
            loop {
                // label 80
                if u0 <= out.u2 {
                    out.ncase = 1;
                    if firstpt && u0 != out.u2 {
                        let mut rff = out.rffu0;
                        out.iu0 = elliptic::ellcubicreal(
                            &p, -out.u1, one, out.u2, -one, out.u3, -one, 0.0, 0.0, &mut rff, u0,
                            out.u2,
                        );
                        out.rffu0 = rff;
                    } else if u0 == out.u2 {
                        out.iu0 = 0.0;
                    }
                    // Table 2 Row 1
                    let m1 = (out.u3 - out.u2) / (out.u3 - out.u1);
                    let jarg =
                        ((out.u3 - out.u1) * dd).sqrt() * half * (imu - su * out.iu0 / dd.sqrt());
                    let (_sn, cn, dn) = sncndn(jarg, m1);
                    let cd2 = cn * cn / (dn * dn);
                    uf = out.u1 + (out.u2 - out.u1) * cd2;
                    let arg = su * (imu - su * out.iu0 / dd.sqrt());
                    tpr = ((fsign(arg) + 1.0) / 2.0) as i32;
                    if pht {
                        let mut rff = out.rffu1;
                        let _ = elliptic::ellcubicreal(
                            &p, -out.u1, one, out.u2, -one, out.u3, -one, 0.0, 0.0, &mut rff, uf,
                            out.u2,
                        );
                        out.rffu1 = rff;
                    }
                    break;
                } else if u0 >= out.u3 {
                    out.ncase = 2;
                    if firstpt && u0 != out.u3 {
                        let mut rff = out.rffu0;
                        out.iu0 = elliptic::ellcubicreal(
                            &p, -out.u1, one, -out.u2, one, -out.u3, one, 0.0, 0.0, &mut rff,
                            out.u3, u0,
                        );
                        out.rffu0 = rff;
                    } else if u0 != out.u3 {
                        // upstream quirk (line 2808): zero the cached value
                        out.iu0 = 0.0;
                    }
                    // Table 2 Row 2
                    let m1 = (out.u3 - out.u2) / (out.u3 - out.u1);
                    let jarg =
                        ((out.u3 - out.u1) * dd).sqrt() * half * (imu + su * out.iu0 / dd.sqrt());
                    let (_sn, cn, dn) = sncndn(jarg, m1);
                    let dc2 = dn * dn / (cn * cn);
                    uf = out.u1 + (out.u3 - out.u1) * dc2;
                    let arg = su * (imu + out.iu0 / dd.sqrt());
                    tpr = ((-fsign(arg) + 1.0) / 2.0) as i32;
                    if pht {
                        let mut rff = out.rffu1;
                        let _ = elliptic::ellcubicreal(
                            &p, -out.u1, one, -out.u2, one, -out.u3, one, 0.0, 0.0, &mut rff,
                            out.u3, uf,
                        );
                        out.rffu1 = rff;
                    }
                    break;
                } else {
                    // unphysical region: modify u0 and retry
                    if su == 1.0 {
                        u0 = out.u3;
                    } else {
                        u0 = out.u2;
                    }
                }
            }
        } else if dis.abs() < 1e-16 {
            // cubic with equal roots
            tpr = 0;
            out.u1 = -two * qq.sqrt() - cc / dd / 3.0;
            out.u2 = -two * qq.sqrt() * (two * pi / 3.0).cos() - cc / dd / 3.0;
            out.u3 = out.u2;
            tpr = 0;
            if u0 <= out.u2 {
                out.ncase = 1;
                if uf > out.u2 {
                    uf = out.u2;
                    // iu = su*1e300 (local, discarded by upstream)
                } else {
                    let sarg = ((u0 - out.u1) / (out.u2 - out.u1)).sqrt();
                    let jarg = ((out.u2 - out.u1) * dd).sqrt() * half * su * imu
                        + half * ((one + sarg) / (one - sarg)).ln();
                    uf = out.u1 + (out.u2 - out.u1) * jarg.tanh().powi(2);
                }
            } else if u0 >= out.u2 {
                out.ncase = 2;
                if uf < out.u2 {
                    uf = out.u2;
                } else {
                    let sarg = ((out.u2 - out.u1) / (u0 - out.u1)).sqrt();
                    let jarg = -((out.u2 - out.u1) * dd).sqrt() * half * su * imu
                        + half * ((one + sarg) / (one - sarg)).ln();
                    uf = out.u1 + (out.u2 - out.u1) / jarg.tanh().powi(2);
                }
            }
        } else {
            // cubic complex case with one real root
            out.ncase = 3;
            tpr = 0;
            let aa = -fsign(rr) * (rr.abs() + dis.sqrt()).powf(third);
            let bb = if aa != 0.0 { qq / aa } else { 0.0 };
            out.u1 = (aa + bb) - cc / dd / 3.0;
            let f = -one / dd / out.u1;
            let g = f / out.u1;
            let mut dummy = 0.0f64;
            let iut = if su > 0.0 {
                elliptic::ellcubiccomplex(
                    &p, -out.u1, one, 0.0, 0.0, f, g, one, &mut dummy, u0, uplus,
                ) / dd.sqrt()
            } else {
                elliptic::ellcubiccomplex(
                    &p, -out.u1, one, 0.0, 0.0, f, g, one, &mut dummy, 0.0, u0,
                ) / dd.sqrt()
            };
            if imu > iut {
                out.uf = -1.0;
                out.ncase = 0;
                out.tpr = 0;
                out.iu = imu;
                return out;
            }
            let _m = -g / 2.0;
            if firstpt {
                out.iu0 =
                    su * elliptic::ellcubiccomplex(
                        &p, -out.u1, one, 0.0, 0.0, f, g, one, &mut dummy, out.u1, u0,
                    ) / dd.sqrt();
            }
            let c3 = if a != 0.0 || l != 0.0 {
                (a + l) / (a - l)
            } else {
                -one
            };
            let c2 = (out.u1 * (three * out.u1 + c3)).sqrt();
            let c1 = (c2 * dd).sqrt();
            let m1 = half + (6.0 * out.u1 + c3) / (8.0 * c2);
            let jarg = c1 * (imu + out.iu0);
            let (_sn, cn, _dn) = sncndn(jarg, m1);
            uf = (c2 + out.u1 - (c2 - out.u1) * cn) / (one + cn);
            if pht {
                let mut rff = out.rffu0;
                let _ = elliptic::ellcubiccomplex(
                    &p, -out.u1, one, 0.0, 0.0, f, g, one, &mut rff, u0, uf,
                );
                out.rffu0 = rff;
            }
        }
    } else if ee == 0.0 && dd == 0.0 {
        // special case q2=0 and l=a: mu=0 at all times
        out.ncase = 4;
        uf = -1.0;
        tpr = 0;
    } else {
        // ----------------------------------------------------------------
        // quartic cases
        // ----------------------------------------------------------------
        // Roots of M(mu) in the biquadratic case.
        let yy = -0.5
            * (a * a - ql2
                + fsign(a * a - ql2) * ((a * a - ql2).powi(2) + 4.0 * q2 * a * a).sqrt());
        let (mut mneg, mut mpos) = if (a * a - ql2) < 0.0 {
            (-yy / a / a, q2 / yy)
        } else {
            (q2 / yy, -yy / a / a)
        };
        let mut muplus = mpos.sqrt();
        let a1 = sm;
        let a2 = sm * neg_one_pow(tpm + 1);
        let a3 = two * f_int((2 * tpm - sm as i32 + 1) as f64 / 4.0) as f64;
        if mneg < 0.0 {
            // symmetric roots case
            if mu0 > muplus {
                muplus = mu0;
                mpos = muplus * muplus;
            }
            if firstpt {
                let vals = imu::calcimusym(a, mneg, mpos, mu0, muplus);
                i1mu = vals.imu0;
                i3mu = vals.imum;
                out.rffmu1 = vals.rf0;
                out.rffmu3 = vals.rfc;
            }
            let vals = imu::calcimusymf(a, mneg, mpos, muf, muplus, i3mu);
            i2mu = vals.imu0;
            out.rffmu2 = vals.rf0;
            imu = a1 * i1mu + a2 * i2mu + a3 * i3mu;
        } else {
            // asymmetric roots case
            if muf.abs() < mneg.sqrt() {
                out.uf = -1.0;
                out.ncase = 0;
                out.tpr = 0;
                out.iu = imu;
                return out;
            } else {
                if fsign(mu0) == -1.0 {
                    muplus = -muplus;
                }
                if muplus.abs() < mu0.abs() {
                    muplus = mu0;
                    mpos = muplus * muplus;
                }
                mneg = (mu0 * mu0).min(mneg);
                if firstpt {
                    let vals = imu::calcimuasym(a, mneg, mpos, mu0, muplus);
                    i1mu = vals.imu0;
                    i3mu = vals.imum;
                    out.rffmu1 = vals.rf0;
                    out.rffmu3 = vals.rfc;
                }
                let vals = imu::calcimuasymf(a, mneg, mpos, muf, muplus, i3mu);
                i2mu = vals.imu0;
                out.rffmu2 = vals.rf0;
            }
            imu = a1 * i1mu + a2 * i2mu + a3 * i3mu;
        }
        // quartic roots
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
                tpr = 0;
                out.u1 = root[0].re;
                if root[1].im == 0.0 {
                    out.u4 = root[1].re;
                } else {
                    out.u4 = root[3].re;
                }
            } else if nreal == 0 {
                out.ncase = 6;
                tpr = 0;
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
                        if su == 1.0 {
                            u0 = out.u3;
                        } else {
                            u0 = out.u2;
                        }
                    }
                }
            }
        }
        if out.ncase == 5 {
            // Table 2 Row 5
            let qs = fsign(q2);
            let f = -qs * one / ee.abs() / out.u1 / out.u4;
            let g = (out.u4 + out.u1) / out.u1 / out.u4 * f;
            let mut dummy = 0.0f64;
            let iut = if su > 0.0 {
                elliptic::ellquarticcomplex(
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
                    &mut dummy,
                    u0,
                    uplus,
                ) / ee.abs().sqrt()
            } else {
                elliptic::ellquarticcomplex(
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
                    &mut dummy,
                    0.0,
                    u0,
                ) / ee.abs().sqrt()
            };
            if imu > iut {
                out.uf = -1.0;
                out.ncase = 0;
                out.tpr = 0;
                out.iu = imu;
                return out;
            }
            let (_ua, ub) = if qs == one {
                (out.u4, out.u1)
            } else {
                (out.u1, out.u4)
            };
            let m = -g * half;
            let n2 = f - g * g / 4.0;
            let c4 = ((m - out.u4).powi(2) + n2).sqrt();
            let c5 = ((m - out.u1).powi(2) + n2).sqrt();
            let c1 = (ee * c4 * c5).abs().sqrt();
            if firstpt {
                out.iu0 =
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
                        &mut dummy,
                        ub,
                        u0,
                    ) / ee.abs().sqrt();
            }
            let m1 = qs * ((c4 + qs * c5).powi(2) - (out.u4 - out.u1).powi(2)) / (4.0 * c4 * c5);
            let jarg = c1 * (imu + out.iu0);
            let (_sn, cn, _dn) = sncndn(jarg, m1);
            uf = (out.u4 * c5 + qs * out.u1 * c4 - (qs * out.u4 * c5 - out.u1 * c4) * cn)
                / ((c4 - qs * c5) * cn + (qs * c4 + c5));
            if pht {
                let mut rff = out.rffu0;
                let _ = elliptic::ellquarticcomplex(
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
                );
                out.rffu0 = rff;
            }
        } else if out.ncase == 6 {
            // Table 2 Row 6: quartic complex with no real roots
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
            let mut h1 = 0.0f64;
            for r in hroots.iter() {
                if r.im == 0.0 {
                    h1 = r.re;
                }
                if h1 != 0.0 {
                    break;
                }
            }
            let h2 = one / h1;
            let g1 = dd / ee / (h2 - h1);
            let g2 = -g1;
            let f1 = one / ee.sqrt();
            let f2 = f1;
            let mut dummy = 0.0f64;
            let iut = if su > 0.0 {
                elliptic::elldoublecomplex(
                    &p, f1, g1, h1, f2, g2, h2, 0.0, 0.0, &mut dummy, u0, uplus,
                ) / ee.abs().sqrt()
            } else {
                elliptic::elldoublecomplex(
                    &p, f1, g1, h1, f2, g2, h2, 0.0, 0.0, &mut dummy, 0.0, u0,
                ) / ee.abs().sqrt()
            };
            if imu > iut {
                out.uf = -1.0;
                out.ncase = 0;
                out.tpr = 0;
                out.iu = imu;
                return out;
            }
            // second root system for mn2
            let coefs2: [Complex64; 7] = [
                Complex64::new(ee.powf(-three), 0.0),
                Complex64::new(-cc / ee.powf(three), 0.0),
                Complex64::new(-ee.powf(-two), 0.0),
                Complex64::new(-ee.powf(-two) * (dd * dd / ee - two * cc), 0.0),
                Complex64::new(-one / ee, 0.0),
                Complex64::new(-cc / ee, 0.0),
                Complex64::new(one, 0.0),
            ];
            let hroots2 = crate::geokerr::special::zroots(&coefs2, 6, true);
            let mut mn2 = 0.0f64;
            for r in hroots2.iter() {
                if r.im == 0.0 {
                    mn2 = r.re;
                }
                if mn2 != 0.0 {
                    break;
                }
            }
            let mut pval = dd / (two * ee.powi(2) * (mn2 * mn2 - one / ee));
            let mut mval = -half * dd / ee - pval;
            if mval < pval {
                std::mem::swap(&mut pval, &mut mval);
            }
            let pr2 = one / (mn2 * ee);
            let n = (mn2 - mval * mval).sqrt();
            let r = (pr2 - pval * pval).sqrt();
            let c4 = ((mval - pval).powi(2) + (n + r).powi(2)).sqrt();
            let c5 = ((mval - pval).powi(2) + (n - r).powi(2)).sqrt();
            let c1 = (c4 + c5) * half * ee.abs().sqrt();
            let c2 = ((4.0 * n * n - (c4 - c5).powi(2)) / ((c4 + c5).powi(2) - 4.0 * n * n)).sqrt();
            let c3 = mval + c2 * n;
            if firstpt {
                out.iu0 = fsign(u0 - c3)
                    * elliptic::elldoublecomplex(
                        &p, f1, g1, h1, f2, g2, h2, 0.0, 0.0, &mut dummy, c3, u0,
                    )
                    / ee.abs().sqrt();
            }
            let m1 = ((c4 - c5) / (c4 + c5)).powi(2);
            let jarg = c1 * (su * imu + out.iu0);
            let (sn, cn, _dn) = sncndn(jarg, m1);
            let sc = sn / cn;
            uf = c3 + (n * (one + c2 * c2) * sc) / (one - c2 * sc);
            if pht {
                let mut rff = out.rffu0;
                let _ = elliptic::elldoublecomplex(
                    &p, f1, g1, h1, f2, g2, h2, 0.0, 0.0, &mut rff, u0, uf,
                );
                out.rffu0 = rff;
            }
        } else {
            // quartic real roots (cases 7, 8)
            if out.ncase == 7 {
                // Table 2 Row 7
                if firstpt && u0 != out.u2 {
                    let mut rff = out.rffu0;
                    out.iu0 = elliptic::ellquarticreal(
                        &p, -out.u1, one, out.u2, -one, out.u3, -one, out.u4, -one, 0.0, 0.0,
                        &mut rff, u0, out.u2,
                    );
                    out.rffu0 = rff;
                }
                let jarg = (ee.abs() * (out.u3 - out.u1) * (out.u4 - out.u2)).sqrt()
                    * half
                    * (imu - su * out.iu0 / (-ee).sqrt());
                let m1 =
                    (out.u4 - out.u1) * (out.u3 - out.u2) / ((out.u4 - out.u2) * (out.u3 - out.u1));
                let (sn, _cn, _dn) = sncndn(jarg, m1);
                let sn2 = sn * sn;
                uf = ((out.u2 - out.u1) * out.u3 * sn2 - out.u2 * (out.u3 - out.u1))
                    / ((out.u2 - out.u1) * sn2 - (out.u3 - out.u1));
                let arg = su * (imu - su * out.iu0 / ee.abs().sqrt());
                tpr = ((fsign(arg) + 1.0) / 2.0) as i32;
                if pht {
                    let mut rff = out.rffu1;
                    let _ = elliptic::ellquarticreal(
                        &p, -out.u1, one, out.u2, -one, out.u3, -one, out.u4, -one, 0.0, 0.0,
                        &mut rff, uf, out.u2,
                    );
                    out.rffu1 = rff;
                }
            } else if out.ncase == 8 {
                // Table 2 Row 8
                if firstpt && u0 != out.u3 {
                    let mut rff = out.rffu0;
                    out.iu0 = elliptic::ellquarticreal(
                        &p, -out.u1, one, -out.u2, one, -out.u3, one, out.u4, -one, 0.0, 0.0,
                        &mut rff, out.u3, u0,
                    );
                    out.rffu0 = rff;
                }
                let jarg = (ee.abs() * (out.u3 - out.u1) * (out.u4 - out.u2)).sqrt()
                    * half
                    * (imu + su * out.iu0 / (-ee).sqrt());
                let m1 =
                    (out.u4 - out.u1) * (out.u3 - out.u2) / ((out.u4 - out.u2) * (out.u3 - out.u1));
                let (sn, _cn, _dn) = sncndn(jarg, m1);
                let sn2 = sn * sn;
                uf = ((out.u4 - out.u3) * out.u2 * sn2 - out.u3 * (out.u4 - out.u2))
                    / ((out.u4 - out.u3) * sn2 - (out.u4 - out.u2));
                let arg = su * (imu + su * out.iu0 / (-ee).sqrt());
                tpr = ((-fsign(arg) + 1.0) / 2.0) as i32;
                if pht {
                    let mut rff = out.rffu1;
                    let _ = elliptic::ellquarticreal(
                        &p, -out.u1, one, -out.u2, one, -out.u3, one, out.u4, -one, 0.0, 0.0,
                        &mut rff, out.u3, uf,
                    );
                    out.rffu1 = rff;
                }
            }
        }
    }

    out.iu = imu;
    out.uf = uf;
    out.tpr = tpr;
    out.i1mu = i1mu;
    out.i2mu = i2mu;
    out.i3mu = i3mu;
    out
}
