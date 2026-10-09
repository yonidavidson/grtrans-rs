//! Carlson elliptic integrals and the cubic/quartic radial integrals.
//!
//! Translation target: `geokerr_wrapper.f` routines RF, RC, RD, RJ,
//! ELLCUBICREAL, ELLCUBICCOMPLEX, ASECH, ELLQUARTICREAL, ELLQUARTICCOMPLEX,
//! ELLDOUBLECOMPLEX, PHIFNKERR, TFNKERR.

use num_complex::Complex64;

/// Carlson symmetric elliptic integral RF(x, y, z).
///
/// Upstream `rf`: R_F = 1/2 \int_0^\infty dt (t+x)^(-1/2)(t+y)^(-1/2)(t+z)^(-1/2).
/// Press et al. (1992), algorithm due to B. C. Carlson; ERRTOL sets accuracy.
pub fn rf(x: f64, y: f64, z: f64) -> f64 {
    const ERRTOL: f64 = 0.0025;
    const THIRD: f64 = 1.0 / 3.0;
    const A1: f64 = 1.0 / 24.0;
    const C2: f64 = 0.1;
    const C3: f64 = 3.0 / 44.0;
    const C4: f64 = 1.0 / 14.0;

    let mut xt = x;
    let mut yt = y;
    let mut zt = z;
    let mut ave = 0.0f64;
    let mut delx = 0.0f64;
    let mut dely = 0.0f64;
    let mut delz = 0.0f64;
    loop {
        let sqrtx = xt.sqrt();
        let sqrty = yt.sqrt();
        let sqrtz = zt.sqrt();
        let alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = THIRD * (xt + yt + zt);
        if ave == 0.0 {
            delx = 0.0;
            dely = 0.0;
            delz = 0.0;
        } else {
            delx = (ave - xt) / ave;
            dely = (ave - yt) / ave;
            delz = (ave - zt) / ave;
        }
        if delx.abs().max(dely.abs()).max(delz.abs()) <= ERRTOL {
            break;
        }
    }
    let e2 = delx * dely - delz * delz;
    let e3 = delx * dely * delz;
    (1.0 + (A1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / ave.sqrt()
}

/// Carlson degenerate elliptic integral RC(x, y).
///
/// Upstream `rc`: R_C = 1/2 \int_0^\infty dt (t+x)^(-1/2)(t+y)^(-1)
/// (Press et al. 1992; algorithm due to B. C. Carlson).
pub fn rc(x: f64, y: f64) -> f64 {
    const ERRTOL: f64 = 0.0012;
    const THIRD: f64 = 1.0 / 3.0;
    const A1: f64 = 0.3;
    const C2: f64 = 1.0 / 7.0;
    const C3: f64 = 0.375;
    const C4: f64 = 9.0 / 22.0;

    let mut xt;
    let mut yt;
    let w;
    if y > 0.0 {
        xt = x;
        yt = y;
        w = 1.0;
    } else {
        xt = x - y;
        yt = -y;
        w = x.sqrt() / xt.sqrt();
    }
    let mut ave;
    let mut s;
    loop {
        let alamb = 2.0 * xt.sqrt() * yt.sqrt() + yt;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        ave = THIRD * (xt + yt + yt);
        s = (yt - ave) / ave;
        if s.abs() <= ERRTOL {
            break;
        }
    }
    w * (1.0 + s * s * (A1 + s * (C2 + s * (C3 + s * C4)))) / ave.sqrt()
}

/// Carlson symmetric elliptic integral RD(x, y, z).
///
/// Upstream `rd`: R_D = 3/2 \int_0^\infty dt (t+x)^(-1/2)(t+y)^(-1/2)(t+z)^(-3/2)
/// (Press et al. 1992; algorithm due to B. C. Carlson).
pub fn rd(x: f64, y: f64, z: f64) -> f64 {
    const ERRTOL: f64 = 0.0015;
    const A1: f64 = 3.0 / 14.0;
    const C2: f64 = 1.0 / 6.0;
    const C3: f64 = 9.0 / 22.0;
    const C4: f64 = 3.0 / 26.0;
    const C5: f64 = 0.25 * C3;
    const C6: f64 = 1.5 * C4;

    let mut xt = x;
    let mut yt = y;
    let mut zt = z;
    let mut sum = 0.0f64;
    let mut fac = 1.0f64;
    let mut ave = 0.0f64;
    let mut delx = 0.0f64;
    let mut dely = 0.0f64;
    let mut delz = 0.0f64;
    loop {
        let sqrtx = xt.sqrt();
        let sqrty = yt.sqrt();
        let sqrtz = zt.sqrt();
        let alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        sum += fac / (sqrtz * (zt + alamb));
        fac *= 0.25;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = 0.2 * (xt + yt + 3.0 * zt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
        if delx.abs().max(dely.abs()).max(delz.abs()) <= ERRTOL {
            break;
        }
    }
    let ea = delx * dely;
    let eb = delz * delz;
    let ec = ea - eb;
    let ed = ea - 6.0 * eb;
    let ee = ed + ec + ec;
    3.0 * sum
        + fac
            * (1.0
                + ed * (-A1 + C5 * ed - C6 * delz * ee)
                + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea)))
            / (ave * ave.sqrt())
}

/// Carlson symmetric elliptic integral RJ(x, y, z, p).
///
/// Upstream `rj`: R_J = 3/2 \int_0^\infty dt
/// (t+x)^(-1/2)(t+y)^(-1/2)(t+z)^(-1/2)(t+p)^(-1)
/// (Press et al. 1992; algorithm due to B. C. Carlson; calls RC and RF).
pub fn rj(x: f64, y: f64, z: f64, p: f64) -> f64 {
    const ERRTOL: f64 = 0.0015;
    const A1: f64 = 3.0 / 14.0;
    const C2: f64 = 1.0 / 3.0;
    const C3: f64 = 3.0 / 22.0;
    const C4: f64 = 3.0 / 26.0;
    const C5: f64 = 0.75 * C3;
    const C6: f64 = 1.5 * C4;
    const C7: f64 = 0.5 * C2;
    const C8: f64 = C3 + C3;

    let mut sum = 0.0f64;
    let mut fac = 1.0f64;
    let mut a = 0.0f64;
    let mut b = 0.0f64;
    let mut rcx = 0.0f64;
    let mut xt;
    let mut yt;
    let mut zt;
    let mut pt;
    if p > 0.0 {
        xt = x;
        yt = y;
        zt = z;
        pt = p;
    } else {
        xt = x.min(y).min(z);
        zt = x.max(y).max(z);
        yt = x + y + z - xt - zt;
        a = 1.0 / (yt - p);
        b = a * (zt - yt) * (yt - xt);
        pt = yt + b;
        let rho = xt * zt / yt;
        let tau = p * pt / yt;
        rcx = rc(rho, tau);
    }
    let mut ave = 0.0f64;
    let mut delx = 0.0f64;
    let mut dely = 0.0f64;
    let mut delz = 0.0f64;
    let mut delp = 0.0f64;
    loop {
        let sqrtx = xt.sqrt();
        let sqrty = yt.sqrt();
        let sqrtz = zt.sqrt();
        let alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        let alpha_t = pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz;
        let alpha = alpha_t * alpha_t;
        let beta_t = pt + alamb;
        let beta = pt * (beta_t * beta_t);
        sum += fac * rc(alpha, beta);
        fac *= 0.25;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        pt = 0.25 * (pt + alamb);
        ave = 0.2 * (xt + yt + zt + pt + pt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
        delp = (ave - pt) / ave;
        if delx.abs().max(dely.abs()).max(delz.abs()).max(delp.abs()) <= ERRTOL {
            break;
        }
    }
    let ea = delx * (dely + delz) + dely * delz;
    let eb = delx * dely * delz;
    let ec = delp * delp;
    let ed = ea - 3.0 * ec;
    let ee = eb + 2.0 * delp * (ea - ec);
    let mut rj = 3.0 * sum
        + fac
            * (1.0
                + ed * (-A1 + C5 * ed - C6 * ee)
                + eb * (C7 + delp * (-C8 + delp * C4))
                + delp * ea * (C2 - delp * C3)
                - C2 * delp * ec)
            / (ave * ave.sqrt());
    if p <= 0.0 {
        rj = a * (b * rj + 3.0 * (rcx - rf(xt, yt, zt)));
    }
    rj
}

/// Integral for Table 1 Rows 1, 2 (upstream `ellcubicreal`).
///
/// Computes \int_y^x dt \Pi_{i=1}^4 (a_i+b_i t)^{p_i/2}. `rff` is the RF
/// piece: written when `p[3] == 0`, read as an input otherwise.
#[allow(clippy::too_many_arguments)]
pub fn ellcubicreal(
    p: &[i32; 5],
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
    a3: f64,
    b3: f64,
    a4: f64,
    b4: f64,
    rff: &mut f64,
    y: f64,
    x: f64,
) -> f64 {
    let one = 1.0f64;
    let half = 0.5f64;
    let two = 2.0f64;
    let three = 3.0f64;

    let mut ellcubic = 0.0f64;
    // (2.1) Carlson (1989)
    let d12 = a1 * b2 - a2 * b1;
    let d13 = a1 * b3 - a3 * b1;
    let d14 = a1 * b4 - a4 * b1;
    let d24 = a2 * b4 - a4 * b2;
    let d34 = a3 * b4 - a4 * b3;
    // (2.2) Carlson (1989)
    let x1 = (a1 + b1 * x).sqrt();
    let x2 = (a2 + b2 * x).sqrt();
    let x3 = (a3 + b3 * x).sqrt();
    let x4 = (a4 + b4 * x).sqrt();
    let y1 = (a1 + b1 * y).sqrt();
    let y2 = (a2 + b2 * y).sqrt();
    let y3 = (a3 + b3 * y).sqrt();
    let y4 = (a4 + b4 * y).sqrt();
    // (2.3) Carlson (1989)
    let u1c = (x1 * y2 * y3 + y1 * x2 * x3) / (x - y);
    let u12 = u1c * u1c;
    let t22 = (x2 * y1 * y3 + y2 * x1 * x3) / (x - y);
    let u22 = t22 * t22;
    let t32 = (x3 * y1 * y2 + y3 * x1 * x2) / (x - y);
    let u32 = t32 * t32;
    // (2.4) Carlson (1989) (upstream first assigns W22=U12, then overwrites it)
    let w22 = u12 - b4 * d12 * d13 / d14;
    // (2.5) Carlson (1989)
    let x4y4 = x4 * y4 / x1 / y1;
    let q22 = x4y4 * x4y4 * w22;
    let p22 = q22 + b4 * d24 * d34 / d14;
    // Now, compute the three integrals we need [-1,-1,-1], [-1,-1,-1,-2], and
    // [-1,-1,-1,-4]:
    if p[3] == 0 {
        // (2.21) Carlson (1989)
        *rff = rf(u32, u22, u12);
        ellcubic = two * *rff;
    } else {
        // (2.12) Carlson (1989)
        let i1c = two * *rff;
        if p[3] == -2 {
            // (2.14) Carlson (1989)
            let i3c =
                two * rc(p22, q22) - two * b1 * d12 * d13 / three / d14 * rj(u32, u22, u12, w22);
            // (2.49) Carlson (1989)
            ellcubic = (b4 * i3c - b1 * i1c) / d14;
        } else {
            // (2.1) Carlson (1989)
            let r12 = a1 / b1 - a2 / b2;
            let r13 = a1 / b1 - a3 / b3;
            let r24i = b2 * b4 / (a2 * b4 - a4 * b2);
            let r34i = b3 * b4 / (a3 * b4 - a4 * b3);
            // (2.13) Carlson (1989)
            let i2c = two / three * d12 * d13 * rd(u32, u22, u12) + two * x1 * y1 / u1c;
            // (2.59) & (2.6) Carlson (1989)
            let k2c =
                b2 * b3 * i2c - two * b4 * (x1 * x2 * x3 / (x4 * x4) - y1 * y2 * y3 / (y4 * y4));
            // (2.62) Carlson (1989)
            let b1d14 = b1 / d14;
            ellcubic = half * b4 / d14 / d24 / d34 * k2c
                + b1d14 * b1d14 * (one - half * r12 * r13 * r24i * r34i) * i1c;
        }
    }
    ellcubic
}

/// Integral for Table 1 Row 3 (upstream `ellcubiccomplex`).
///
/// Computes \int_y^x dt \Pi_{i=1,4} (a_i+b_i t)^{p_i/2} (f+gt+ht^2)^{p_2/2}.
/// `rff` is the RF piece: written when `p[3] == 0`, read as an input otherwise.
#[allow(clippy::too_many_arguments)]
pub fn ellcubiccomplex(
    p: &[i32; 5],
    a1: f64,
    b1: f64,
    a4: f64,
    b4: f64,
    f: f64,
    g: f64,
    h: f64,
    rff: &mut f64,
    y: f64,
    x: f64,
) -> f64 {
    let one = 1.0f64;
    let two = 2.0f64;
    let half = 0.5f64;
    let three = 3.0f64;
    let four = 4.0f64;
    let six = 6.0f64;

    let mut ellcubic = 0.0f64;
    let x1 = (a1 + b1 * x).sqrt();
    let x4 = (a4 + b4 * x).sqrt();
    let y1 = (a1 + b1 * y).sqrt();
    let y4 = (a4 + b4 * y).sqrt();
    let d14 = a1 * b4 - a4 * b1;
    // (2.2) Carlson (1991)
    let beta1 = g * b1 - two * h * a1;
    let _beta4 = g * b4 - two * h * a4; // computed upstream, unused
                                        // (2.3) Carlson (1991)
    let a11 = (two * f * b1 * b1 - two * g * a1 * b1 + two * h * a1 * a1).sqrt();
    let c44 = (two * f * b4 * b4 - two * g * a4 * b4 + two * h * a4 * a4).sqrt();
    let a142 = two * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + two * h * a1 * a4;
    // (2.4) Carlson (1991)
    let xi = (f + g * x + h * x * x).sqrt();
    let eta = (f + g * y + h * y * y).sqrt();
    // (3.1) Carlson (1991)
    let m2t = (x1 + y1) * ((xi + eta) * (xi + eta) - h * (x - y) * (x - y)).sqrt() / (x - y);
    let m2 = m2t * m2t;
    // (3.2) Carlson (1991)
    let s2h = (two * h).sqrt();
    let lp2 = m2 - beta1 + s2h * a11;
    let lm2 = m2 - beta1 - s2h * a11;
    if p[3] == 0 {
        *rff = rf(m2, lm2, lp2);
        // (1.2) Carlson (1991)
        ellcubic = four * *rff;
    } else {
        // (3.8) Carlson (1991)
        let i1c = four * *rff;
        // (3.3) Carlson (1991)
        let u = (x1 * eta + y1 * xi) / (x - y);
        let u2 = u * u;
        let wp2 = m2 - b1 * (a142 + a11 * c44) / d14;
        let a11sq = a11 * a11;
        let w2 = u2 - a11sq * b4 / two / d14;
        // (3.4) Carlson (1991)
        let x4y4 = x4 * y4 / x1 / y1;
        let q2 = x4y4 * x4y4 * w2;
        let c44sq = c44 * c44;
        let p2 = q2 + c44sq * b4 / two / d14;
        // (3.5) Carlson (1991)
        let rho = s2h * a11 - beta1;
        // (3.9) Carlson (1991)
        if p[3] == -2 {
            // (2.49) Carlson (1989)
            let i3c = (two * a11 / three / c44)
                * ((-four * b1 / d14) * (a142 + a11 * c44) * rj(m2, lm2, lp2, wp2) - six * *rff
                    + three * rc(u2, w2))
                + two * rc(p2, q2);
            ellcubic = (b4 * i3c - b1 * i1c) / d14;
        } else {
            // (2.19) Carlson (1991)
            let r24xr34 = half * c44sq / h / (b4 * b4);
            let r12xr13 = half * a11sq / h / (b1 * b1);
            // (3.11) Carlson (1991)
            let n2c = two / three * s2h / a11
                * (four * rho * rd(m2, lm2, lp2) - six * *rff + three / u)
                + two / x1 / y1 / u;
            // (2.5) & (3.12) Carlson (1991)
            let k2c = half * a11sq * n2c - two * d14 * (xi / x1 / (x4 * x4) - eta / y1 / (y4 * y4));
            // (2.62) Carlson (1989)
            let b1d14 = b1 / d14;
            ellcubic = half / d14 / (h * b4 * r24xr34) * k2c
                + b1d14 * b1d14 * (one - half * r12xr13 / r24xr34) * i1c;
        }
    }
    ellcubic
}

/// Inverse hyperbolic secant (upstream `asech`).
///
/// Upstream returns 0 for x <= 0.
pub fn asech(x: f64) -> f64 {
    if x > 0.0 {
        ((1.0 + (1.0 - x * x).sqrt()) / x).ln()
    } else {
        0.0
    }
}

/// Integral for Table 1 Rows 7, 8 (upstream `ellquarticreal`).
///
/// Computes \int_y^x dt \Pi_{i=1}^5 (a_i+b_i t)^{p_i/2}. `rff` is the RF
/// piece: written when `p[4] == 0`, read as an input otherwise.
#[allow(clippy::too_many_arguments)]
pub fn ellquarticreal(
    p: &[i32; 5],
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
    a3: f64,
    b3: f64,
    a4: f64,
    b4: f64,
    a5: f64,
    b5: f64,
    rff: &mut f64,
    y: f64,
    x: f64,
) -> f64 {
    let one = 1.0f64;
    let half = 0.5f64;
    let two = 2.0f64;
    let three = 3.0f64;

    // (2.1) Carlson (1988)
    let d12 = a1 * b2 - a2 * b1;
    let d13 = a1 * b3 - a3 * b1;
    let d14 = a1 * b4 - a4 * b1;
    let d24 = a2 * b4 - a4 * b2;
    let d34 = a3 * b4 - a4 * b3;
    let d15 = a1 * b5 - a5 * b1;
    let d25 = a2 * b5 - a5 * b2;
    let d35 = a3 * b5 - a5 * b3;
    let d45 = a4 * b5 - a5 * b4;
    // (2.2) Carlson (1988)
    let x1 = (a1 + b1 * x).sqrt();
    let x2 = (a2 + b2 * x).sqrt();
    let x3 = (a3 + b3 * x).sqrt();
    let x4 = (a4 + b4 * x).sqrt();
    let y1 = (a1 + b1 * y).sqrt();
    let y2 = (a2 + b2 * y).sqrt();
    let y3 = (a3 + b3 * y).sqrt();
    let y4 = (a4 + b4 * y).sqrt();
    // (2.3) Carlson (1988)
    let u122t = (x1 * x2 * y3 * y4 + y1 * y2 * x3 * x4) / (x - y);
    let u122 = u122t * u122t;
    let u132t = (x1 * x3 * y2 * y4 + y1 * y3 * x2 * x4) / (x - y);
    let u132 = u132t * u132t;
    let u142t = (x1 * x4 * y2 * y3 + y1 * y4 * x2 * x3) / (x - y);
    let u142 = u142t * u142t;
    // Now, compute the three integrals we need [-1,-1,-1,-1],[-1,-1,-1,-1,-2],
    // [-1,-1,-1,-1,-4]:
    if p[4] == 0 {
        *rff = rf(u122, u132, u142);
        // (2.17) Carlson (1988)
        two * *rff
    } else {
        // (2.13) Carlson (1988)
        let i1 = two * *rff;
        let x52 = a5 + b5 * x;
        let y52 = a5 + b5 * y;
        // (2.4) Carlson (1988)
        let w22 = u122 - d13 * d14 * d25 / d15;
        // (2.5) Carlson (1989)
        let x1y1 = x1 * y1;
        let q22 = x52 * y52 / (x1y1 * x1y1) * w22;
        let p22 = q22 + d25 * d35 * d45 / d15;
        // (2.15) Carlson (1988)
        if p[4] == -2 {
            // (2.35) Carlson (1988)
            let i3 = two * d12 * d13 * d14 / three / d15 * rj(u122, u132, u142, w22)
                + two * rc(p22, q22);
            (b5 * i3 - b1 * i1) / d15
        } else {
            let i2 = two / three * d12 * d13 * rd(u122, u132, u142)
                + two * x1 * y1 / x4 / y4 / u142.sqrt();
            // (2.1) Carlson (1988)
            let r12 = a1 / b1 - a2 / b2;
            let r13 = a1 / b1 - a3 / b3;
            let r25i = b2 * b5 / (a2 * b5 - a5 * b2);
            let r35i = b3 * b5 / (a3 * b5 - a5 * b3);
            // (2.48) Carlson (1988)
            let a111m1m2 = x1 * x2 * x3 / x4 / x52 - y1 * y2 * y3 / y4 / y52;
            let b5sq = b5 * b5;
            let b1d15 = b1 / d15;
            half * b5sq * d24 * d34 / d15 / d25 / d35 / d45 * i2
                + b1d15 * b1d15 * (one - half * r12 * r13 * r25i * r35i) * i1
                - b5sq / d15 / d25 / d35 * a111m1m2
        }
        // The upstream trailing `ellquarticreal = 0; return` is unreachable.
    }
}

/// Integral for Table 1 Row 5 (upstream `ellquarticcomplex`).
///
/// Computes \int_y^x dt \Pi_{i=1,4,5} (a_i+b_i t)^{p_i/2} (f+gt+ht^2)^{p_2/2}.
/// `rff` is the RF piece: written when `p[4] == 0`, read as an input otherwise.
#[allow(clippy::too_many_arguments)]
pub fn ellquarticcomplex(
    p: &[i32; 5],
    a1: f64,
    b1: f64,
    a4: f64,
    b4: f64,
    a5: f64,
    b5: f64,
    f: f64,
    g: f64,
    h: f64,
    rff: &mut f64,
    y: f64,
    x: f64,
) -> f64 {
    let one = 1.0f64;
    let two = 2.0f64;
    let half = 0.5f64;
    let three = 3.0f64;
    let four = 4.0f64;
    let six = 6.0f64;

    // (2.1) Carlson (1991)
    let x1 = (a1 + b1 * x).sqrt();
    let x4 = (a4 + b4 * x).sqrt();
    let y1 = (a1 + b1 * y).sqrt();
    let y4 = (a4 + b4 * y).sqrt();
    // (2.3) Carlson (1991)
    let a11 = (two * f * b1 * b1 - two * g * a1 * b1 + two * h * a1 * a1).sqrt();
    let c44 = (two * f * b4 * b4 - two * g * a4 * b4 + two * h * a4 * a4).sqrt();
    let a142 = two * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + two * h * a1 * a4;
    // (2.4) Carlson (1991)
    let xi = (f + g * x + h * x * x).sqrt();
    let eta = (f + g * y + h * y * y).sqrt();
    // (2.6) Carlson (1991)
    let m2t =
        (x1 * y4 + y1 * x4) * ((xi + eta) * (xi + eta) - h * (x - y) * (x - y)).sqrt() / (x - y);
    let m2 = m2t * m2t;
    // (2.7) Carlson (1991)
    let lp2 = m2 + a142 + a11 * c44;
    let lm2 = (m2 + a142 - a11 * c44).max(0.0);
    if p[4] == 0 {
        // (1.2) Carlson (1991)
        *rff = rf(m2, lm2, lp2);
        four * *rff
    } else {
        // (2.14) Carlson (1991)
        let i1 = four * *rff;
        // (2.1) Carlson (1991)
        let d14 = a1 * b4 - a4 * b1;
        let d15 = a1 * b5 - a5 * b1;
        let d45 = a4 * b5 - a5 * b4;
        let x52 = a5 + b5 * x;
        let y52 = a5 + b5 * y;
        // (2.3) Carlson (1991)
        let c552 = two * f * b5 * b5 - two * g * a5 * b5 + two * h * a5 * a5;
        let c55 = c552.sqrt();
        let a152 = two * f * b1 * b5 - g * (a1 * b5 + a5 * b1) + two * h * a5 * a1;
        // (2.8) Carlson (1991)
        let u = (x1 * x4 * eta + y1 * y4 * xi) / (x - y);
        let u2 = u * u;
        let wp2 = m2 + d14 * (a152 + a11 * c55) / d15;
        let a11sq = a11 * a11;
        let w2 = u2 - a11sq * d45 / two / d15;
        // (2.9) Carlson (1991)
        let x1y1 = x1 * y1;
        let q2 = x52 * y52 / (x1y1 * x1y1) * w2;
        let c55sq = c55 * c55;
        let p2 = q2 + c55sq * d45 / two / d15;
        // (2.15) Carlson (1991)
        if p[4] == -2 {
            // (2.35) Carlson (1988)
            let i3 = (two * a11 / three / c55)
                * ((four * d14 / d15) * (a152 + a11 * c55) * rj(m2, lm2, lp2, wp2) - six * *rff
                    + three * rc(u2, w2))
                + two * rc(p2, q2);
            (b5 * i3 - b1 * i1) / d15
        } else {
            let i2 = two * a11 / three / c44
                * (four * (a142 + a11 * c44) * rd(m2, lm2, lp2) - six * *rff)
                + two * a11 / c44 / u
                + two * x1 * y1 / x4 / y4 / u;
            // (2.5) Carlson (1991)
            let a111m1m2 = x1 * xi / x4 / x52 - y1 * eta / y4 / y52;
            let b5sq = b5 * b5;
            let b1sq = b1 * b1;
            let c44sq = c44 * c44;
            b5sq / (two * d15 * d45) * c44sq / c552 * i2
                + b1sq / (d15 * d15) * (one - half * b5sq * a11sq / b1sq / c552) * i1
                - two * b5sq / d15 / c552 * a111m1m2
        }
    }
}

/// Integral for Table 1 Row 6 (upstream `elldoublecomplex`).
///
/// Computes \int_y^x dt (f_1+g_1t+h_1t^2)^{p_1/2} (f_2+g_2t+h_2t^2)^{p_2/2}
/// (a_5+b_5t)^{p_5/2} following Carlson (1992). `rff` is the RF piece:
/// written when `p[4] == 0`, read as an input otherwise.
#[allow(clippy::too_many_arguments)]
pub fn elldoublecomplex(
    p: &[i32; 5],
    f1: f64,
    g1: f64,
    h1: f64,
    f2: f64,
    g2: f64,
    h2: f64,
    a5: f64,
    b5: f64,
    rff: &mut f64,
    y: f64,
    x: f64,
) -> f64 {
    let one = 1.0f64;
    let half = 0.5f64;
    let two = 2.0f64;
    let three = 3.0f64;
    let four = 4.0f64;
    let _six = 6.0f64;

    // (2.1) Carlson (1992)
    let xi1 = (f1 + g1 * x + h1 * x * x).sqrt();
    let xi2 = (f2 + g2 * x + h2 * x * x).sqrt();
    let eta1 = (f1 + g1 * y + h1 * y * y).sqrt();
    let eta2 = (f2 + g2 * y + h2 * y * y).sqrt();
    // (2.4) Carlson (1992)
    let theta1 = two * f1 + g1 * (x + y) + two * h1 * x * y;
    let theta2 = two * f2 + g2 * (x + y) + two * h2 * x * y;
    // (2.5) Carlson (1992)
    let zeta1 = (2.0 * xi1 * eta1 + theta1).sqrt();
    let zeta2 = (2.0 * xi2 * eta2 + theta2).sqrt();
    // (2.6) Carlson (1992)
    let m = zeta1 * zeta2 / (x - y);
    let m2 = m * m;
    // (2.7) Carlson (1992)
    let delta122 = two * f1 * h2 + two * f2 * h1 - g1 * g2;
    let delta112 = four * f1 * h1 - g1 * g1;
    let delta222 = four * f2 * h2 - g2 * g2;
    let delta = (delta122 * delta122 - delta112 * delta222).sqrt();
    // (2.8) Carlson (1992)
    let deltap = delta122 + delta;
    let deltam = delta122 - delta;
    let lp2 = m2 + deltap;
    let lm2 = m2 + deltam;
    if p[4] == 0 {
        // (2.36) Carlson (1992)
        *rff = rf(m2, lm2, lp2);
        four * *rff
    } else {
        // (2.6) Carlson (1992)
        let u = (xi1 * eta2 + eta1 * xi2) / (x - y);
        let u2 = u * u;
        // (2.11) Carlson (1992)
        let alpha15 = two * f1 * b5 - g1 * a5;
        let alpha25 = two * f2 * b5 - g2 * a5;
        let beta15 = g1 * b5 - two * h1 * a5;
        let beta25 = g2 * b5 - two * h2 * a5;
        // (2.12) Carlson (1992)
        let gamma1 = half * (alpha15 * b5 - beta15 * a5);
        let gamma2 = half * (alpha25 * b5 - beta25 * a5);
        // (2.13) Carlson (1992)
        let lambda = delta112 * gamma2 / gamma1;
        let omega2 = m2 + lambda;
        let psi = half * (alpha15 * beta25 - alpha25 * beta15);
        let psi2 = psi * psi;
        // (2.15) Carlson (1992)
        let xi5 = a5 + b5 * x;
        let eta5 = a5 + b5 * y;
        // (2.16) Carlson (1992)
        let am111m1 = one / xi1 * xi2 - one / eta1 * eta2;
        let a1111m4 = xi1 * xi2 / (xi5 * xi5) - eta1 * eta2 / (eta5 * eta5);
        // (2.17) Carlson (1992)
        let xx =
            xi5 * eta5 * (theta1 * half * am111m1 - xi5 * eta5 * a1111m4) / ((x - y) * (x - y));
        // (2.18) Carlson (1992)
        let s = half * (m2 + delta122) - u2;
        let s2 = s * s;
        // (2.19) Carlson (1992)
        let mu = gamma1 * xi5 * eta5 / xi1 / eta1;
        let t = mu * s + two * gamma1 * gamma2;
        let t2 = t * t;
        let v2 = mu * mu * (s2 + lambda * u2);
        // (2.20) Carlson (1992)
        let b2 = omega2 * omega2 * (s2 / u2 + lambda);
        let a2 = b2 + lambda * lambda * psi2 / gamma1 / gamma2;
        // (2.22) Carlson (1992)
        let h = delta112 * psi * (rj(m2, lm2, lp2, omega2) / three + half * rc(a2, b2))
            / (gamma1 * gamma1)
            - xx * rc(t2, v2);
        if p[4] == -2 {
            // (2.39) Carlson (1992)
            -two * (b5 * h + beta15 * *rff / gamma1)
        } else {
            let a1111m2 = xi1 * xi2 / xi5 - eta1 * eta2 / eta5;
            // (2.2) Carlson (1992)
            let xi1p = half * (g1 + two * h1 * x) / xi1;
            let eta1p = half * (g1 + two * h1 * y) / eta1;
            // (2.3) Carlson (1992)
            let b = xi1p * xi2 - eta1p * eta2;
            // (2.9) Carlson (1992)
            let g = two / three * delta * deltap * rd(m2, lm2, lp2)
                + half * delta / u
                + (delta122 * theta1 - delta112 * theta2) / four / xi1 / eta1 / u;
            // (2.10) Carlson (1992)
            let sigma = g - deltap * *rff + b;
            // (2.41) Carlson (1992)
            b5 * (beta15 / gamma1 + beta25 / gamma2) * h
                + beta15 * beta15 * *rff / (gamma1 * gamma1)
                + b5 * b5 * (sigma - b5 * a1111m2) / gamma1 / gamma2
        }
    }
}

/// Azimuthal integrand contribution (upstream `phifnkerr`).
///
/// Closed-form antiderivative of the azimuthal equation for the equal-roots
/// cubic case; `u1 < u2 = u3`.
pub fn phifnkerr(u: f64, u1: f64, u2: f64, l: f64, a: f64) -> f64 {
    let two = 2.0f64;
    let up = 1.0 / (1.0 + (1.0 - a * a).sqrt());
    let umi = 1.0 - (1.0 - a * a).sqrt();
    let s = (u - u1).sqrt();
    let su21 = (u2 - u1).sqrt();
    let _su1 = (-u1).sqrt(); // computed upstream, unused
    let sup = (1.0 - u1 / up).sqrt();
    let sumv = (1.0 - u1 * umi).sqrt();

    (l + two * a * u2 - two * l * u2) / (su21 * (u2 * umi - 1.0) * (u2 / up - 1.0))
        * ((su21 + s) / (su21 - s)).ln()
        + (two * a + l * (-two + umi))
            * umi.sqrt()
            * ((sumv + s * umi.sqrt()) / (sumv - s * umi.sqrt())).ln()
            / (sumv * (u2 * umi - 1.0) * (umi - 1.0 / up))
        + (two * a + l * (-two + 1.0 / up))
            * (1.0 / up).sqrt()
            * ((sup + (1.0 / up).sqrt() * s) / (sup - (1.0 / up).sqrt() * s)).ln()
            / (1.0 / up - umi)
            * sup
            * (u2 / up - 1.0)
}

/// Time integrand contribution (upstream `tfnkerr`).
///
/// Closed-form antiderivative of the time equation for the equal-roots cubic
/// case; `u1 < u2 = u3`.
pub fn tfnkerr(u0: f64, uf: f64, u1: f64, u2: f64, l: f64, a: f64) -> f64 {
    let two = 2.0f64;
    let up = 1.0 / (1.0 + (1.0 - a * a).sqrt());
    let umi = 1.0 - (1.0 - a * a).sqrt();
    let sf = (uf - u1).sqrt();
    let ss = (u0 - u1).sqrt();
    let su21 = (u2 - u1).sqrt();
    let su1 = (-u1).sqrt();
    let sup = (1.0 - u1 / up).sqrt();
    let sumv = (1.0 - u1 * umi).sqrt();
    let u2sq = u2 * u2;
    let u2cu = u2sq * u2;
    let umicu = umi * umi * umi;

    sf / (uf * u1 * u2)
        - ss / (u0 * u1 * u2)
        - ((-u2 - two * u1 * (1.0 + u2 * umi + u2 / up))
            * (((sf + su1) * (su1 - ss)) / ((-sf + su1) * (ss + su1))).ln())
            / (two * (-u1).powf(3.0 / 2.0) * u2sq)
        + ((1.0 - two * a * l * u2cu + a * a * u2sq * (1.0 + two * u2))
            * (((su21 + sf) * (su21 - ss)) / ((su21 - sf) * (su21 + ss))).ln())
            / (u2sq * su21 * (u2 * umi - 1.0) * (u2 / up - 1.0))
        + ((-two * a * l + a * a * (two + umi) + umicu)
            * umi.sqrt()
            * (((sumv + sf * umi.sqrt()) * (sumv - ss * umi.sqrt()))
                / ((sumv - sf * umi.sqrt()) * (sumv + ss * umi.sqrt())))
            .ln())
            / (sumv * (u2 * umi - 1.0) * (umi - 1.0 / up))
        + ((-two * a * l + a * a * (two + 1.0 / up) + 1.0 / (up * up * up))
            * (1.0 / up).sqrt()
            * (((sup + sf * (1.0 / up).sqrt()) * (sup - ss * (1.0 / up).sqrt()))
                / ((sup - sf * (1.0 / up).sqrt()) * (sup + ss * (1.0 / up).sqrt())))
            .ln())
            / ((1.0 / up - umi) * sup * (u2 / up - 1.0))
}

/// Re-export for callers that need complex roots.
pub type C64 = Complex64;

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(got: f64, want: f64, tol: f64, what: &str) {
        let scale = want.abs().max(1.0);
        assert!(
            (got - want).abs() <= tol * scale,
            "{what}: got {got}, want {want}"
        );
    }

    #[test]
    fn rf_matches_analytic_values() {
        for &x in &[0.5f64, 1.0, 2.0, 10.0] {
            assert_close(rf(x, x, x), x.powf(-0.5), 1e-12, "rf(x,x,x)");
        }
    }

    #[test]
    fn rd_matches_analytic_values() {
        for &x in &[0.5f64, 1.0, 2.0, 10.0] {
            assert_close(rd(x, x, x), x.powf(-1.5), 1e-12, "rd(x,x,x)");
        }
    }

    #[test]
    fn rj_matches_analytic_values() {
        for &x in &[0.5f64, 1.0, 2.0, 10.0] {
            assert_close(rj(x, x, x, x), x.powf(-1.5), 1e-12, "rj(x,x,x,x)");
        }
    }

    #[test]
    fn rc_matches_analytic_values() {
        assert_close(rc(1.0, 1.0), 1.0, 1e-14, "rc(1,1)");
        for &x in &[0.5f64, 1.0, 2.0, 10.0] {
            assert_close(rc(x, x), x.powf(-0.5), 1e-12, "rc(x,x)");
        }
    }

    #[test]
    fn asech_matches_analytic_values() {
        let want = ((1.0 + (0.75f64).sqrt()) / 0.5).ln();
        assert_close(asech(0.5), want, 1e-15, "asech(0.5)");
        assert_eq!(asech(0.0), 0.0);
        assert_eq!(asech(-0.5), 0.0);
        assert_eq!(asech(1.0), 0.0);
    }

    #[test]
    fn ellcubicreal_elementary_case() {
        // p(4)=0, a_i=b_i=1:
        //   \int_y^x dt / (1+t)^{3/2} = 2 (1/sqrt(1+y) - 1/sqrt(1+x)).
        let p = [-1i32, -1, -1, 0, 0];
        let (y, x) = (2.0f64, 4.0f64);
        let mut rff = 0.0;
        let got = ellcubicreal(&p, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &mut rff, y, x);
        let want = 2.0 * (1.0 / (1.0 + y).sqrt() - 1.0 / (1.0 + x).sqrt());
        assert_close(got, want, 1e-13, "ellcubicreal");
        assert_close(rff, 0.5 * want, 1e-13, "ellcubicreal rff");
    }

    #[test]
    fn ellquarticreal_elementary_case() {
        // p(5)=0, a_i=b_i=1:
        //   \int_y^x dt / (1+t)^2 = 1/(1+y) - 1/(1+x).
        let p = [-1i32, -1, -1, -1, 0];
        let (y, x) = (2.0f64, 4.0f64);
        let mut rff = 0.0;
        let got = ellquarticreal(
            &p, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &mut rff, y, x,
        );
        let want = 1.0 / (1.0 + y) - 1.0 / (1.0 + x);
        assert_close(got, want, 1e-13, "ellquarticreal");
        assert_close(rff, 0.5 * want, 1e-13, "ellquarticreal rff");
    }
}
