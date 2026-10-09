//! Coordinate time and azimuth along the geodesic.
//!
//! Translation target: `geokerr_wrapper.f` routines GEOPHITIME (lines
//! 1648-2357), CALCPHITMUASYM (2361-2432), CALCPHITMUSYM (2435-2504) and
//! ELLPHITMU (2508-2590).
//!
//! References to "Eq. (x)" are to Dexter & Agol (2009).

use crate::geokerr::elliptic::{
    ellcubiccomplex, ellcubicreal, elldoublecomplex, ellquarticcomplex, ellquarticreal, phifnkerr,
    rd, rj, tfnkerr,
};

/// State carried between GEOPHITIME calls (computed once on the first point
/// of a geodesic: RDC, RJC, TU01..TU04, TMU1, TMU3, PHIMU1, PHIMU3).
#[derive(Clone, Copy, Debug, Default)]
pub struct PhitimeState {
    pub rdc: f64,
    pub rjc: f64,
    pub tu01: f64,
    pub tu02: f64,
    pub tu03: f64,
    pub tu04: f64,
    pub tmu1: f64,
    pub tmu3: f64,
    pub phimu1: f64,
    pub phimu3: f64,
}

/// Outputs of one GEOPHITIME call.
#[derive(Clone, Copy, Debug, Default)]
pub struct PhitimeOutput {
    pub phimu: f64,
    pub tmu: f64,
    pub phiu: f64,
    pub tu: f64,
    pub lambda: f64,
    pub state: PhitimeState,
}

/// Fortran `SIGN(1.d0, x)`: +1 for `x >= 0`, -1 otherwise. Negative zero
/// compares equal to zero, matching the Fortran `SIGN` definition (this
/// differs from `f64::copysign` on `-0.0`).
#[inline]
fn fsgn1(x: f64) -> f64 {
    if x >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// Fortran `(-1)**n` for an integer exponent, evaluated by parity (valid for
/// negative `n` too, though upstream only uses non-negative counters).
#[inline]
fn neg1_pow(n: i32) -> f64 {
    if n % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

/// Outputs of `ellphitmu` (upstream's output dummy arguments). `rdc`/`rjc`
/// are also inputs when `firstpt == false`.
#[derive(Clone, Copy, Debug, Default)]
struct EllphitmuOut {
    elle0: f64,
    ellef: f64,
    ellec: f64,
    ellpi0: f64,
    ellpif: f64,
    ellpic: f64,
    rdc: f64,
    rjc: f64,
}

/// Upstream `ellphitmu` (lines 2508-2590): Legendre elliptic integrals of the
/// 1st/2nd/3rd kind for the MU0, MUF and MUPLUS pieces. `n` is in the
/// Numerical Recipes convention; `rdc`/`rjc` are inputs when
/// `firstpt == false`.
#[allow(clippy::too_many_arguments)]
fn ellphitmu(
    phi0: f64,
    phif: f64,
    m1: f64,
    n: f64,
    firstpt: bool,
    rf0: f64,
    rff: f64,
    rfc: f64,
    rdc_in: f64,
    rjc_in: f64,
) -> EllphitmuOut {
    let pi = (-1.0_f64).acos();
    let mut rdc = rdc_in;
    let mut rjc = rjc_in;
    let mut elle0 = 0.0;
    let mut ellec = 0.0;
    let mut ellpi0 = 0.0;
    let mut ellpic = 0.0;
    // Convert to argument k from complimentary parameter m1 (line 2543)
    let ak = (1.0 - m1).sqrt();
    // First do mu0, muplus components if they haven't already been calculated
    if firstpt {
        let s0 = phi0.sin();
        let q0 = (1.0 - s0 * ak) * (1.0 + s0 * ak);
        // Carlson integrals for Legendre's of the 2nd kind using Press et al
        // (1992) 6.11.20 (lines 2549-2552)
        let rd0 = rd(1.0 - s0 * s0, q0, 1.0);
        rdc = rd(0.0, m1, 1.0);
        elle0 = s0 * (rf0 - ((s0 * ak) * (s0 * ak)) * rd0 / 3.0);
        ellec = 2.0 * (rfc - ak * ak * rdc / 3.0);
        // n only equals 1 if it was set that way, in which case the phi terms
        // are zero (lines 2554-2563)
        if n != 1.0 {
            // Legendre's 3rd kind from Press et al (1992) 6.11.21
            let rj0 = rj(1.0 - s0 * s0, q0, 1.0, 1.0 - n * s0 * s0);
            rjc = rj(0.0, m1, 1.0, 1.0 - n);
            ellpi0 = s0 * (rf0 + n * s0 * s0 * rj0 / 3.0);
            ellpic = 2.0 * (rfc + n * rjc / 3.0);
        } else {
            ellpi0 = 0.0;
            ellpic = 0.0;
        }
        // If phi > pi/2, we've computed int_0^pi/2 - int_pi/2^phi; use
        // int_0^phi = 2 int_0^pi/2 - (int_0^pi/2 - int_pi/2^phi) (lines 2566-2569)
        if phi0 > pi / 2.0 {
            ellpi0 = ellpic - ellpi0;
            elle0 = ellec - elle0;
        }
    }
    // Repeat procedure for muf part (lines 2572-2587)
    let sf = phif.sin();
    let qf = (1.0 - sf * ak) * (1.0 + sf * ak);
    let rdf = rd(1.0 - sf * sf, qf, 1.0);
    let mut ellef = sf * (rff - ((sf * ak) * (sf * ak)) * rdf / 3.0);
    let mut ellpif = 0.0;
    if n != 1.0 {
        let rjf = rj(1.0 - sf * sf, qf, 1.0, 1.0 - n * sf * sf);
        ellpif = sf * (rff + n * sf * sf * rjf / 3.0);
    } else {
        ellpif = 0.0;
    }
    // See above note on phi > pi/2.
    if phif > pi / 2.0 {
        ellpic = 2.0 * (rfc + n * rjc / 3.0);
        ellec = 2.0 * (rfc - ak * ak * rdc / 3.0);
        ellpif = ellpic - ellpif;
        ellef = ellec - ellef;
    }
    EllphitmuOut {
        elle0,
        ellef,
        ellec,
        ellpi0,
        ellpif,
        ellpic,
        rdc,
        rjc,
    }
}

/// Outputs of `calcphitmuasym`/`calcphitmusym`. `phimu0`/`tmu0` etc. are only
/// computed when `firstpt == true`; otherwise they are pass-throughs.
#[derive(Clone, Copy, Debug, Default)]
struct CalcphitmuOut {
    phimu0: f64,
    phimuf: f64,
    phimum: f64,
    tmu0: f64,
    tmuf: f64,
    tmum: f64,
    rdc: f64,
    rjc: f64,
}

/// Upstream `calcphitmuasym` (lines 2361-2432): integral pieces of phimu and
/// tmu in the asymmetric (M_- > 0) case.
#[allow(clippy::too_many_arguments)]
fn calcphitmuasym(
    a: f64,
    mneg: f64,
    mpos: f64,
    mu0: f64,
    muf: f64,
    muplus: f64,
    rf0: f64,
    rff: f64,
    rfc: f64,
    rdc_in: f64,
    rjc_in: f64,
    firstpt: bool,
) -> CalcphitmuOut {
    // After Eq. (45) (lines 2400-2401)
    let fac = a.abs() * muplus;
    let m1 = mneg / mpos;
    // If the orbit crosses the pole, don't compute phi integral terms
    // (lines 2403-2412)
    let (n, fac2) = if mpos != 1.0 {
        // After Eq. (46)
        ((mpos - mneg) / (1.0 - mpos), 1.0 - mpos)
    } else {
        // Set n to a value we can detect
        (-1.0, 1.0)
    };
    let phi0 = ((mpos - mu0 * mu0) / (mpos - mneg)).sqrt().asin();
    let phif = ((mpos - muf * muf) / (mpos - mneg)).sqrt().asin();
    // Get Legendre integrals of 1st/3rd kind (line 2416)
    let out = ellphitmu(phi0, phif, m1, -n, firstpt, rf0, rff, rfc, rdc_in, rjc_in);
    let mut res = CalcphitmuOut {
        rdc: out.rdc,
        rjc: out.rjc,
        ..Default::default()
    };
    // Now calculate required integrals (lines 2418-2430)
    if firstpt {
        // Eq. (45) for the mu0 term and the integral between the turning
        // points. The factors of two in tmum and phimum set the integral
        // between turning points equal to the complete integral.
        res.tmu0 = fac * out.elle0;
        res.tmum = fac * out.ellec / 2.0;
        // Eq. (46) for the mu0 term and the integral between turning points.
        res.phimu0 = out.ellpi0 / fac2 / fac;
        res.phimum = out.ellpic / fac2 / fac / 2.0;
    }
    // Eq. (45) for the muf term
    res.tmuf = fac * out.ellef;
    // Eq. (46) for the muf term.
    res.phimuf = out.ellpif / fac2 / fac;
    res
}

/// Upstream `calcphitmusym` (lines 2435-2504): integral pieces of phimu and
/// tmu in the symmetric (M_- < 0) case.
#[allow(clippy::too_many_arguments)]
fn calcphitmusym(
    a: f64,
    mneg: f64,
    mpos: f64,
    mu0: f64,
    muf: f64,
    muplus: f64,
    rf0: f64,
    rff: f64,
    rfc: f64,
    rdc_in: f64,
    rjc_in: f64,
    firstpt: bool,
) -> CalcphitmuOut {
    // Immediately after Eq. (45) (lines 2472-2473)
    let fac = a.abs() * (mpos - mneg).sqrt();
    let m1 = -mneg / (mpos - mneg);
    // If the orbit crosses the pole, don't compute phi integral terms
    // (lines 2475-2484)
    let (n, fac2) = if mpos != 1.0 {
        // Immediately after Eq. (46)
        (mpos / (1.0 - mpos), 1.0 - mpos)
    } else {
        // Set n to a value we can detect
        (-1.0, 1.0)
    };
    // phi arguments for Legendre integrals (in the paper, x = sin(phi))
    let phi0 = (mu0 / muplus).acos();
    let phif = (muf / muplus).acos();
    // Use ellpi to calculate other pieces if non-zero (line 2489)
    let out = ellphitmu(phi0, phif, m1, -n, firstpt, rf0, rff, rfc, rdc_in, rjc_in);
    let mut res = CalcphitmuOut {
        rdc: out.rdc,
        rjc: out.rjc,
        ..Default::default()
    };
    // Calculate integrals involving mu0 and turning points only if this is
    // the first point on the geodesic (lines 2491-2498)
    if firstpt {
        // Eq. (45) for the mu0 term and the integral between turning points.
        res.tmu0 = fac * out.elle0;
        res.tmum = fac * out.ellec;
        // Eq. (46) for the mu0 term and the integral between turning points.
        res.phimu0 = out.ellpi0 / fac2 / fac;
        res.phimum = out.ellpic / fac2 / fac;
    }
    // Eq. (45) for the muf term.
    res.tmuf = fac * out.ellef;
    // Eq. (46) for the muf term.
    res.phimuf = out.ellpif / fac2 / fac;
    res
}

/// Upstream `geophitime`: compute delta phi and delta t between
/// (U0, MU0) and (UF, MUF), plus the affine parameter. `state_in` carries
/// the first-point data (only used when `firstpt == false`).
#[allow(clippy::too_many_arguments)]
pub fn geophitime(
    u0: f64,
    uf: f64,
    mu0: f64,
    muf: f64,
    a: f64,
    l: f64,
    l2: f64,
    q2: f64,
    tpm: i32,
    tpr: i32,
    su: f64,
    sm: f64,
    iu: f64,
    h1: f64,
    u1: f64,
    u2: f64,
    u3: f64,
    u4: f64,
    rffu0: f64,
    rffu1: f64,
    rffmu1: f64,
    rffmu2: f64,
    rffmu3: f64,
    ncase: i32,
    state_in: &PhitimeState,
    firstpt: bool,
) -> PhitimeOutput {
    // line 1714: pi = acos(-1)
    let pi = (-1.0_f64).acos();
    // lines 1716-1718
    let mut p = [-1_i32, -1, -1, 0, 0];
    // line 1719
    let uplus = 1.0 / (1.0 + (1.0 - a * a).sqrt());
    // line 1721: only calculate 1/u_- since u_- blows up when a=0
    let umi = 1.0 - (1.0 - a * a).sqrt();
    // line 1722
    let ur = -1.0 / (2.0 * (1.0 - a * a).sqrt());
    // line 1723
    let qs = fsgn1(q2);
    // line 1724, Eq. (22)
    let dd = 2.0 * ((a - l) * (a - l) + q2);
    // line 1726, Eq. (22)
    let ee = -a * a * q2;
    // line 1727
    let ql2 = q2 + l2;
    // line 1728
    let mut a1 = sm;
    // line 1729
    let mut a2 = sm * neg1_pow(tpm);
    // line 1731, Eq. (35)
    let mut a3 = 2.0 * ((2.0 * f64::from(tpm) + 3.0 - sm) / 4.0).trunc() - 1.0;

    // State variables. These mirror upstream's in/out dummy arguments: they
    // are computed when firstpt=true and otherwise start from state_in, but
    // upstream can write them on later calls as well, so they are always
    // written back into the returned state.
    let mut rdc = state_in.rdc;
    let mut rjc = state_in.rjc;
    let mut tu01 = state_in.tu01;
    let mut tu02 = state_in.tu02;
    let mut tu03 = state_in.tu03;
    let mut tu04 = state_in.tu04;
    let mut tmu1 = state_in.tmu1;
    let mut tmu3 = state_in.tmu3;
    let mut phimu1 = state_in.phimu1;
    let mut phimu3 = state_in.phimu3;

    // Return values and local temporaries (upstream declares several locals
    // without initializing them; we start them at zero to stay deterministic).
    let mut tu = 0.0_f64;
    let mut phiu = 0.0_f64;
    let mut phimu = 0.0_f64;
    let mut tmu = 0.0_f64;
    let mut lambdau = 0.0_f64;
    let mut tmu2 = 0.0_f64;
    let mut phimu2 = 0.0_f64;
    let mut tu0 = 0.0_f64;
    let mut tu1 = 0.0_f64;
    let mut tu2 = 0.0_f64;
    let mut tu3 = 0.0_f64;
    let mut tu4 = 0.0_f64;
    let mut phiu0 = 0.0_f64;
    let mut phiu1 = 0.0_f64;
    let mut phiu2 = 0.0_f64;
    let mut tu11 = 0.0_f64;
    let mut tu12 = 0.0_f64;
    let mut tu13 = 0.0_f64;
    let mut tu14 = 0.0_f64;
    let mut phiu01 = 0.0_f64;
    let mut phiu02 = 0.0_f64;
    let mut phiu11 = 0.0_f64;
    let mut phiu12 = 0.0_f64;
    // RF pieces are in/out arguments upstream; keep them mutable so any
    // write-back (p=0 calls) would be visible to later calls as in Fortran.
    let mut rffu0 = rffu0;
    let mut rffu1 = rffu1;
    let mut muplus;

    if ncase == 0 {
        // line 1732-1737
        tu = 0.0;
        phiu = 0.0;
        phimu = 0.0;
        tmu = 0.0;
        lambdau = 0.0;
    } else if ncase < 3 {
        // These are the cubic real roots cases with u1<0<u2<=u3 (line 1738)
        p[4] = 0;
        if (u3 - u2).abs() < 1.0e-12 {
            // These are the equal roots cases (lines 1741-1748)
            tu = su * tfnkerr(u0, uf, u1, u2, l, a) / dd.sqrt();
            phiu = su * (phifnkerr(uf, u1, u2, l, a) - phifnkerr(u0, u1, u2, l, a)) / dd.sqrt();
            if u0 >= u3 {
                phiu = -phiu;
                tu = -tu;
            }
        } else if u0 <= u2 {
            // Table 1 Row 1 (line 1749)
            if firstpt && u0 != u2 {
                p[3] = -2;
                // First three u0 integrals in Eq. (47)
                phiu01 = ellcubicreal(
                    &p,
                    -u1,
                    1.0,
                    u2,
                    -1.0,
                    u3,
                    -1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu0,
                    u0,
                    u2,
                );
                phiu02 = ellcubicreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, 1.0, -umi, &mut rffu0, u0, u2,
                );
                tu02 = ellcubicreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, 0.0, 1.0, &mut rffu0, u0, u2,
                );
                tu03 = phiu01;
                tu04 = phiu02;
                p[3] = -4;
                // Fourth u0 integral in Eq. (47)
                tu01 = ellcubicreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, 0.0, 1.0, &mut rffu0, u0, u2,
                );
            } else if u0 == u2 {
                tu01 = 0.0;
                tu02 = 0.0;
                tu03 = 0.0;
                tu04 = 0.0;
                phiu01 = 0.0;
                phiu02 = 0.0;
            } else {
                phiu01 = tu03;
                phiu02 = tu04;
            }
            if (uf - u2).abs() > 1.0e-16 {
                p[3] = -2;
                // First three uf integrals in Eq. (47)
                phiu11 = ellcubicreal(
                    &p,
                    -u1,
                    1.0,
                    u2,
                    -1.0,
                    u3,
                    -1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu1,
                    uf,
                    u2,
                );
                phiu12 = ellcubicreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, 1.0, -umi, &mut rffu1, uf, u2,
                );
                tu12 = ellcubicreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, 0.0, 1.0, &mut rffu1, uf, u2,
                );
                tu13 = phiu11;
                tu14 = phiu12;
                p[3] = -4;
                // Fourth uf integral in Eq. (47)
                tu11 = ellcubicreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, 0.0, 1.0, &mut rffu1, uf, u2,
                );
            } else {
                tu11 = 0.0;
                tu12 = 0.0;
                tu13 = 0.0;
                tu14 = 0.0;
                phiu11 = 0.0;
                phiu12 = 0.0;
            }
            // Eq. (47) for u0 piece
            tu0 = (umi - 1.0 / uplus) * tu01 + (umi * umi - 1.0 / (uplus * uplus)) * tu02
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu03
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu04;
            // Eq. (47) for uf piece
            tu1 = (umi - 1.0 / uplus) * tu11 + (umi * umi - 1.0 / (uplus * uplus)) * tu12
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu13
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu14;
            // Combine according to Eq. (54)
            tu = su * (tu0 - neg1_pow(tpr) * tu1) / dd.sqrt() * ur;
            // Eq. (48) for u0 piece
            phiu0 = -(l / uplus + 2.0 * (a - l)) * phiu01 + (l * umi + 2.0 * (a - l)) * phiu02;
            // Eq. (48) for uf piece
            phiu1 = -(l / uplus + 2.0 * (a - l)) * phiu11 + (l * umi + 2.0 * (a - l)) * phiu12;
            // Combine according to Eq. (54)
            phiu = su * (phiu0 - neg1_pow(tpr) * phiu1) / dd.sqrt() * ur;
            // Eq. (49)
            lambdau = su * (tu01 - neg1_pow(tpr) * tu11) / dd.sqrt();
        } else if u0 >= u3 {
            // Table 1 Row 2 (line 1810)
            if firstpt && u0 != u3 {
                p[3] = -2;
                // First three u0 integrals in Eq. (47)
                phiu01 = -ellcubicreal(
                    &p,
                    -u1,
                    1.0,
                    -u2,
                    1.0,
                    -u3,
                    1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu0,
                    u3,
                    u0,
                );
                phiu02 = -ellcubicreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, 1.0, -umi, &mut rffu0, u3, u0,
                );
                tu02 = -ellcubicreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, 0.0, 1.0, &mut rffu0, u3, u0,
                );
                tu03 = phiu01;
                tu04 = phiu02;
                // Fourth u0 integral in Eq. (47)
                p[3] = -4;
                tu01 = -ellcubicreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, 0.0, 1.0, &mut rffu0, u3, u0,
                );
            } else if u0 == u3 {
                tu01 = 0.0;
                tu02 = 0.0;
                tu03 = 0.0;
                tu04 = 0.0;
                phiu01 = 0.0;
                phiu02 = 0.0;
            } else {
                phiu01 = tu03;
                phiu02 = tu04;
            }
            if uf != u3 {
                p[3] = -2;
                // First three uf integrals in Eq. (47)
                phiu11 = -ellcubicreal(
                    &p,
                    -u1,
                    1.0,
                    -u2,
                    1.0,
                    -u3,
                    1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu1,
                    u3,
                    uf,
                );
                phiu12 = -ellcubicreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, 1.0, -umi, &mut rffu1, u3, uf,
                );
                tu12 = -ellcubicreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, 0.0, 1.0, &mut rffu1, u3, uf,
                );
                tu13 = phiu11;
                tu14 = phiu12;
                p[3] = -4;
                // Fourth uf integral in Eq. (47)
                tu11 = -ellcubicreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, 0.0, 1.0, &mut rffu1, u3, uf,
                );
            } else {
                tu11 = 0.0;
                tu12 = 0.0;
                tu13 = 0.0;
                tu14 = 0.0;
                phiu11 = 0.0;
                phiu12 = 0.0;
            }
            // Eq. (47) for u0 part
            tu0 = (umi - 1.0 / uplus) * tu01 + (umi * umi - 1.0 / (uplus * uplus)) * tu02
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu03
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu04;
            // Eq. (47) for uf part
            tu1 = (umi - 1.0 / uplus) * tu11 + (umi * umi - 1.0 / (uplus * uplus)) * tu12
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu13
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu14;
            // Combine according to Eq. (54)
            tu = su * (tu0 - neg1_pow(tpr) * tu1) / dd.sqrt() * ur;
            // Eq. (48) for u0 part
            phiu0 = -(l / uplus + 2.0 * (a - l)) * phiu01 + (l * umi + 2.0 * (a - l)) * phiu02;
            // Eq. (48) for uf part
            phiu1 = -(l / uplus + 2.0 * (a - l)) * phiu11 + (l * umi + 2.0 * (a - l)) * phiu12;
            // Combine according to Eq. (54)
            phiu = su * (phiu0 - neg1_pow(tpr) * phiu1) / dd.sqrt() * ur;
            // Eq. (49)
            lambdau = su * (tu01 - neg1_pow(tpr) * tu11) / dd.sqrt();
        }
    } else if ncase == 3 {
        // This is a cubic complex case with one real root. Table 1 Row 3
        // (line 1872)
        let f = -1.0 / dd / u1;
        let g = f / u1;
        let h = 1.0;
        if u0 < uf {
            p[3] = -4;
            // Fourth integral in Eq. (47)
            tu1 = ellcubiccomplex(&p, -u1, 1.0, 0.0, 1.0, f, g, h, &mut rffu0, u0, uf);
            p[3] = -2;
            // First three integrals in Eq. (47)
            tu2 = ellcubiccomplex(&p, -u1, 1.0, 0.0, 1.0, f, g, h, &mut rffu0, u0, uf);
            tu3 = ellcubiccomplex(&p, -u1, 1.0, 1.0, -1.0 / uplus, f, g, h, &mut rffu0, u0, uf);
            tu4 = ellcubiccomplex(&p, -u1, 1.0, 1.0, -umi, f, g, h, &mut rffu0, u0, uf);
            // Eq. (47)
            tu = su * ur / dd.sqrt()
                * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                    - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu3
                    + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
            // Minus sign in phi components is from flipping the sign of the
            // u/u_\pm-1 factor to keep arguments positive.
            phiu0 = -tu3;
            phiu1 = -tu4;
            // Eq. (48)
            phiu = su * ur / dd.sqrt()
                * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1);
            lambdau = su * tu1 / dd.sqrt();
        } else if u0 > uf {
            p[3] = -4;
            // Fourth integral in Eq. (47)
            tu1 = -ellcubiccomplex(&p, -u1, 1.0, 0.0, 1.0, f, g, h, &mut rffu0, uf, u0);
            p[3] = -2;
            // First three integrals in Eq. (47)
            tu2 = -ellcubiccomplex(&p, -u1, 1.0, 0.0, 1.0, f, g, h, &mut rffu0, uf, u0);
            tu3 = -ellcubiccomplex(&p, -u1, 1.0, 1.0, -1.0 / uplus, f, g, h, &mut rffu0, uf, u0);
            tu4 = -ellcubiccomplex(&p, -u1, 1.0, 1.0, -umi, f, g, h, &mut rffu0, uf, u0);
            // Eq. (47)
            tu = su * ur / dd.sqrt()
                * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                    - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu3
                    + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
            // Minus sign in phi components is from flipping the sign of the
            // u/u_\pm-1 factor to keep arguments positive.
            phiu0 = -tu3;
            phiu1 = -tu4;
            // Eq. (48)
            phiu = su * ur / dd.sqrt()
                * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1);
            // Eq. (49)
            lambdau = su * tu1 / dd.sqrt();
        } else {
            tu = 0.0;
            phiu = 0.0;
            lambdau = 0.0;
        }
    }
    if q2 == 0.0 {
        // Calculate mu components in the special case q2=0, where the t and
        // phi integrals are elementary (lines 1924-1940)
        let s1 = fsgn1(mu0);
        a1 = s1 * sm;
        a2 = s1 * sm * neg1_pow(tpm + 1);
        if l.abs() < a.abs() {
            muplus = s1 * (1.0 - l2 / a / a).sqrt();
            phimu1 = fsgn1(muplus * l / a)
                * (((muplus - ((mu0 / muplus).asin() / 2.0).tan())
                    / (1.0 - muplus * muplus).sqrt())
                .atan()
                    + ((muplus + ((mu0 / muplus).asin() / 2.0).tan())
                        / (1.0 - muplus * muplus).sqrt())
                    .atan());
            phimu2 = fsgn1(muplus * l / a)
                * (((muplus - ((muf / muplus).asin() / 2.0).tan())
                    / (1.0 - muplus * muplus).sqrt())
                .atan()
                    + ((muplus + ((muf / muplus).asin() / 2.0).tan())
                        / (1.0 - muplus * muplus).sqrt())
                    .atan());
            phimu = a1 * phimu1 + a2 * phimu2;
            tmu = (muplus * a).abs()
                * (a2 * (1.0 - muf * muf / (muplus * muplus)).sqrt()
                    + a1 * (1.0 - mu0 * mu0 / (muplus * muplus)).sqrt());
        } else {
            phimu = 0.0;
            tmu = 0.0;
        }
    } else if a == 0.0 {
        // Calculate mu components in the special case a=0. The t and phi
        // integrals are again elementary (lines 1941-1965)
        tmu = 0.0;
        muplus = (q2 / ql2).sqrt();
        let mut vfm = (muf - muplus * muplus) / (1.0 - muf) / muplus;
        let mut vfp = -(muf + muplus * muplus) / (1.0 + muf) / muplus;
        // This is code to suppress floating point errors in phimu calculation.
        if vfm.abs() > 1.0 {
            vfm = fsgn1(vfm);
        }
        if vfp.abs() > 1.0 {
            vfp = fsgn1(vfp);
        }
        let mut vsm = if (mu0 - muplus * muplus) == 0.0 {
            -1.0
        } else {
            (mu0 - muplus * muplus) / (1.0 - mu0) / muplus
        };
        let mut vsp = if (mu0 + muplus * muplus) == 0.0 {
            -1.0
        } else {
            -(mu0 + muplus * muplus) / (1.0 + mu0) / muplus
        };
        if vsp.abs() > 1.0 {
            vsp = fsgn1(vsp);
        }
        if vsm.abs() > 1.0 {
            vsm = fsgn1(vsm);
        }
        phimu1 = pi - vsm.asin() + vsp.asin();
        phimu2 = vfm.asin() - vfp.asin() + pi;
        phimu3 = 2.0 * pi;
        phimu = -l * iu + fsgn1(l) * 0.5 * (a1 * phimu1 + a2 * phimu2 + a3 * phimu3);
    }
    if ncase == 4 {
        // This is the special case where q2=0 and l=a. U(u)=1 and the t and
        // phi u components are elementary (lines 1967-1977)
        // Eq. (48)
        phiu = su
            * ur
            * (l * ((uf / uplus - 1.0) / (u0 / uplus - 1.0)).ln()
                - l * umi * umi * ((uf * umi - 1.0) / (u0 * umi - 1.0)).ln());
        // Eq. (47)
        tu = su
            * ur
            * ((umi - 1.0 / uplus) * (1.0 / u0 - 1.0 / uf)
                + (umi * umi - 1.0 / (uplus * uplus)) * (uf / u0).ln()
                + (a * a / uplus + 1.0 / (uplus * uplus * uplus))
                    * uplus
                    * ((uf / uplus - 1.0) / (u0 / uplus - 1.0)).ln()
                - (a * a * umi + umi * umi * umi)
                    * umi
                    * ((uf * umi - 1.0) / (u0 * umi - 1.0)).ln());
        // Eq. (49)
        lambdau = su * (1.0 / u0 - 1.0 / uf);
    } else if ncase == 5 {
        // This is the quartic case with one pair of complex roots
        // (lines 1978-2029). Table 1 Row 5.
        p[3] = -1;
        let f = -qs * 1.0 / ee.abs() / u1 / u4;
        let g = (u4 + u1) / u1 / u4 * f;
        let h = 1.0;
        let apply48 = |phiu0: f64, phiu1: f64| -> f64 {
            su * ur / (-qs * ee).sqrt()
                * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1)
        };
        if u0 < uf {
            p[4] = -4;
            // Fourth integral in Eq. (47)
            tu1 = ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                0.0,
                1.0,
                f,
                g,
                h,
                &mut rffu0,
                u0,
                uf,
            );
            p[4] = -2;
            // First three integrals in Eq. (47)
            tu2 = ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                0.0,
                1.0,
                f,
                g,
                h,
                &mut rffu0,
                u0,
                uf,
            );
            tu3 = ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                1.0,
                -1.0 / uplus,
                f,
                g,
                h,
                &mut rffu0,
                u0,
                uf,
            );
            tu4 = ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                1.0,
                -umi,
                f,
                g,
                h,
                &mut rffu0,
                u0,
                uf,
            );
            // Eq. (47)
            tu = su * ur / (-qs * ee).sqrt()
                * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                    - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu3
                    + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
            // Eq. (49)
            lambdau = su * tu1 / ee.abs().sqrt();
            // Minus sign in phi components is from flipping the sign of the
            // u/u_\pm-1 factor to keep arguments positive.
            phiu0 = -tu3;
            phiu1 = -tu4;
            // Eq. (48)
            phiu = apply48(phiu0, phiu1);
        } else if u0 > uf {
            p[4] = -4;
            // Fourth integral in Eq. (47)
            tu1 = -ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                0.0,
                1.0,
                f,
                g,
                h,
                &mut rffu0,
                uf,
                u0,
            );
            p[4] = -2;
            // First three integrals in Eq. (47)
            tu2 = -ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                0.0,
                1.0,
                f,
                g,
                h,
                &mut rffu0,
                uf,
                u0,
            );
            tu3 = -ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                1.0,
                -1.0 / uplus,
                f,
                g,
                h,
                &mut rffu0,
                uf,
                u0,
            );
            tu4 = -ellquarticcomplex(
                &p,
                -u1,
                1.0,
                u4 * qs,
                -qs,
                1.0,
                -umi,
                f,
                g,
                h,
                &mut rffu0,
                uf,
                u0,
            );
            // Eq. (47)
            tu = su * ur / (-qs * ee).sqrt()
                * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                    - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu3
                    + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
            // Minus sign in phi components is from flipping the sign of the
            // u/u_\pm-1 factor to keep arguments positive.
            phiu0 = -tu3;
            phiu1 = -tu4;
            // Eq. (48)
            phiu = apply48(phiu0, phiu1);
            // Eq. (49)
            lambdau = su * tu1 / ee.abs().sqrt();
        } else {
            tu = 0.0;
            phiu = 0.0;
            lambdau = 0.0;
        }
    } else if ncase == 6 {
        // This is the quartic complex case with no real roots
        // (lines 2030-2082). Table 1 Row 6.
        p[3] = -1;
        let h2 = 1.0 / h1;
        let g1 = dd / ee / (h2 - h1);
        let g2 = -g1;
        let f1 = 1.0 / ee.sqrt();
        let f2 = f1;
        if u0 < uf {
            p[4] = -4;
            // Fourth integral in Eq. (47)
            tu1 = elldoublecomplex(&p, f1, g1, h1, f2, g2, h2, 0.0, 1.0, &mut rffu0, u0, uf);
            p[4] = -2;
            // First three integrals in Eq. (47)
            tu2 = elldoublecomplex(&p, f1, g1, h1, f2, g2, h2, 0.0, 1.0, &mut rffu0, u0, uf);
            tu3 = elldoublecomplex(
                &p,
                f1,
                g1,
                h1,
                f2,
                g2,
                h2,
                1.0,
                -1.0 / uplus,
                &mut rffu0,
                u0,
                uf,
            );
            tu4 = elldoublecomplex(&p, f1, g1, h1, f2, g2, h2, 1.0, -umi, &mut rffu0, u0, uf);
            // Minus sign in phi components is from flipping the sign of the
            // u/u_\pm-1 factor to keep arguments positive.
            phiu0 = -tu3;
            phiu1 = -tu4;
            // Eq. (47)
            tu = su * ur / (-qs * ee).sqrt()
                * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                    - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu3
                    + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
            // Eq. (48)
            phiu = su * ur / (-qs * ee).sqrt()
                * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1);
            // Eq. (49)
            lambdau = su * tu1 / (-qs * ee).sqrt();
        } else if u0 > uf {
            p[4] = -4;
            // Fourth integral in Eq. (47)
            tu1 = -elldoublecomplex(&p, f1, g1, h1, f2, g2, h2, 0.0, 1.0, &mut rffu0, uf, u0);
            p[4] = -2;
            // First three integrals in Eq. (47)
            tu2 = -elldoublecomplex(&p, f1, g1, h1, f2, g2, h2, 0.0, 1.0, &mut rffu0, uf, u0);
            tu3 = -elldoublecomplex(
                &p,
                f1,
                g1,
                h1,
                f2,
                g2,
                h2,
                1.0,
                -1.0 / uplus,
                &mut rffu0,
                uf,
                u0,
            );
            tu4 = -elldoublecomplex(&p, f1, g1, h1, f2, g2, h2, 1.0, -umi, &mut rffu0, uf, u0);
            // Minus sign in phi components is from flipping the sign of the
            // u/u_\pm-1 factor to keep arguments positive.
            phiu0 = -tu3;
            phiu1 = -tu4;
            // Eq. (47)
            tu = su * ur / (-qs * ee).sqrt()
                * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                    - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu3
                    + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
            // Eq. (48)
            phiu = su * ur / (-qs * ee).sqrt()
                * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1);
            // Eq. (49)
            lambdau = su * tu1 / (-qs * ee).sqrt();
        } else {
            tu = 0.0;
            phiu = 0.0;
        }
    } else if ncase > 6 {
        // These are the quartic cases with all real roots (lines 2083-2296)
        p[3] = -1;
        if (u3 - u2).abs() < 1.0e-12 {
            // These are the equal roots quartic cases
            if u0 < uf {
                if u0 < u2 {
                    // Table 1 Row 7
                    p[4] = -2;
                    // First three integrals in Eq. (47)
                    phiu1 = ellquarticreal(
                        &p,
                        -u1,
                        1.0,
                        u2,
                        -1.0,
                        u3,
                        -1.0,
                        u4,
                        -1.0,
                        1.0,
                        -1.0 / uplus,
                        &mut rffu0,
                        u0,
                        uf,
                    );
                    phiu2 = ellquarticreal(
                        &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 1.0, -umi, &mut rffu0, u0, uf,
                    );
                    tu2 = ellquarticreal(
                        &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u0, uf,
                    );
                    tu3 = phiu1;
                    tu4 = phiu2;
                    p[4] = -4;
                    // Fourth integral in Eq. (47)
                    tu1 = ellquarticreal(
                        &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u0, uf,
                    );
                } else if u0 > u3 {
                    // Table 1 Row 8. Note: upstream does not set p(5) before
                    // the first three calls here (it reads the uninitialized
                    // local), whereas the mirrored u0 > uf branch does.
                    // The literal translation leaves p(5) at its initial
                    // value; this path is only reachable for degenerate
                    // (double) roots inside the forbidden region.
                    // First three integrals in Eq. (47)
                    phiu1 = ellquarticreal(
                        &p,
                        -u1,
                        1.0,
                        -u2,
                        1.0,
                        -u3,
                        1.0,
                        u4,
                        -1.0,
                        1.0,
                        -1.0 / uplus,
                        &mut rffu0,
                        u0,
                        uf,
                    );
                    phiu2 = ellquarticreal(
                        &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 1.0, -umi, &mut rffu0, u0, uf,
                    );
                    tu2 = ellquarticreal(
                        &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u0, uf,
                    );
                    tu3 = phiu1;
                    tu4 = phiu2;
                    p[4] = -4;
                    // Fourth integral in Eq. (47)
                    tu1 = ellquarticreal(
                        &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u0, uf,
                    );
                } else {
                    tu11 = 0.0;
                    tu12 = 0.0;
                    tu13 = 0.0;
                    tu14 = 0.0;
                    phiu11 = 0.0;
                    phiu12 = 0.0;
                }
                // Eq. (47)
                tu = su * ur / (-qs * ee).sqrt()
                    * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                        - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus))
                            * tu3
                        + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
                // Eq. (48). Upstream uses phiu0/phiu1 here although this
                // branch only assigned phiu1/phiu2 (and tu3/tu4); both are
                // uninitialized locals upstream, kept at zero here.
                phiu = su * ur / (-qs * ee).sqrt()
                    * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1);
                // Eq. (49)
                lambdau = tu1 / (-qs * ee).sqrt();
            } else if u0 > uf {
                if u0 < u2 {
                    // Table 1 Row 7
                    p[4] = -2;
                    // First three integrals in Eq. (47)
                    phiu1 = -ellquarticreal(
                        &p,
                        -u1,
                        1.0,
                        u2,
                        -1.0,
                        u3,
                        -1.0,
                        u4,
                        -1.0,
                        1.0,
                        -1.0 / uplus,
                        &mut rffu0,
                        uf,
                        u0,
                    );
                    phiu2 = -ellquarticreal(
                        &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 1.0, -umi, &mut rffu0, uf, u0,
                    );
                    tu2 = -ellquarticreal(
                        &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, uf, u0,
                    );
                    tu3 = -phiu1;
                    tu4 = -phiu2;
                    p[4] = -4;
                    // Fourth integral in Eq. (47)
                    tu1 = -ellquarticreal(
                        &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, uf, u0,
                    );
                } else if u0 > u3 {
                    // Table 1 Row 8
                    p[4] = -2;
                    // First three integrals in Eq. (47)
                    phiu1 = -ellquarticreal(
                        &p,
                        -u1,
                        1.0,
                        -u2,
                        1.0,
                        -u3,
                        1.0,
                        u4,
                        -1.0,
                        1.0,
                        -1.0 / uplus,
                        &mut rffu0,
                        uf,
                        u0,
                    );
                    phiu2 = -ellquarticreal(
                        &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 1.0, -umi, &mut rffu0, uf, u0,
                    );
                    tu2 = -ellquarticreal(
                        &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, uf, u0,
                    );
                    tu3 = -phiu1;
                    tu4 = -phiu2;
                    p[4] = -4;
                    // Fourth integral in Eq. (47)
                    tu1 = -ellquarticreal(
                        &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, uf, u0,
                    );
                } else {
                    tu11 = 0.0;
                    tu12 = 0.0;
                    tu13 = 0.0;
                    tu14 = 0.0;
                    phiu11 = 0.0;
                    phiu12 = 0.0;
                }
                // Eq. (47)
                tu = su * ur / (-qs * ee).sqrt()
                    * ((umi - 1.0 / uplus) * tu1 + (umi * umi - 1.0 / (uplus * uplus)) * tu2
                        - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus))
                            * tu3
                        + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu4);
                // Eq. (48)
                phiu = su * ur / (-qs * ee).sqrt()
                    * ((l / uplus + 2.0 * (a - l)) * phiu0 - (l * umi + 2.0 * (a - l)) * phiu1);
                // Eq. (49)
                lambdau = tu1 / (-qs * ee).sqrt();
            } else {
                tu = 0.0;
                phiu = 0.0;
                lambdau = 0.0;
            }
        } else if u0 <= u2 {
            // This is the quartic case with distinct, real roots.
            // Table 1 Row 7 (line 2174)
            if firstpt && u0 != u2 {
                // Only compute integrals involving u0 once per geodesic
                p[4] = -2;
                // First three u0 integrals in Eq. (47)
                phiu01 = ellquarticreal(
                    &p,
                    -u1,
                    1.0,
                    u2,
                    -1.0,
                    u3,
                    -1.0,
                    u4,
                    -1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu0,
                    u0,
                    u2,
                );
                phiu02 = ellquarticreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 1.0, -umi, &mut rffu0, u0, u2,
                );
                tu02 = ellquarticreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u0, u2,
                );
                tu03 = phiu01;
                tu04 = phiu02;
                p[4] = -4;
                // Fourth u0 integral in Eq. (48)
                tu01 = ellquarticreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u0, u2,
                );
            } else if u0 == u2 {
                tu01 = 0.0;
                tu02 = 0.0;
                tu03 = 0.0;
                tu04 = 0.0;
                phiu01 = 0.0;
                phiu02 = 0.0;
            } else {
                phiu01 = tu03;
                phiu02 = tu04;
            }
            if (uf - u2).abs() > 1.0e-12 {
                p[4] = -2;
                // First three uf integrals in Eq. (47)
                phiu11 = ellquarticreal(
                    &p,
                    -u1,
                    1.0,
                    u2,
                    -1.0,
                    u3,
                    -1.0,
                    u4,
                    -1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu1,
                    uf,
                    u2,
                );
                phiu12 = ellquarticreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 1.0, -umi, &mut rffu1, uf, u2,
                );
                tu12 = ellquarticreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu1, uf, u2,
                );
                tu13 = phiu11;
                tu14 = phiu12;
                p[4] = -4;
                // Fourth uf integral in Eq. (47)
                tu11 = ellquarticreal(
                    &p, -u1, 1.0, u2, -1.0, u3, -1.0, u4, -1.0, 0.0, 1.0, &mut rffu1, uf, u2,
                );
            } else {
                tu11 = 0.0;
                tu12 = 0.0;
                tu13 = 0.0;
                tu14 = 0.0;
                phiu11 = 0.0;
                phiu12 = 0.0;
            }
            // Eq. (47) for u0 piece
            tu0 = (umi - 1.0 / uplus) * tu01 + (umi * umi - 1.0 / (uplus * uplus)) * tu02
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu03
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu04;
            // Eq. (47) for uf piece
            tu1 = (umi - 1.0 / uplus) * tu11 + (umi * umi - 1.0 / (uplus * uplus)) * tu12
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu13
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu14;
            tu = su * (tu0 - neg1_pow(tpr) * tu1) / ee.abs().sqrt() * ur;
            // Eq. (48)
            phiu0 = -(l / uplus + 2.0 * (a - l)) * phiu01 + (l * umi + 2.0 * (a - l)) * phiu02;
            phiu1 = -(l / uplus + 2.0 * (a - l)) * phiu11 + (l * umi + 2.0 * (a - l)) * phiu12;
            // Combine according to Eq. (54)
            phiu = su * (phiu0 - neg1_pow(tpr) * phiu1) / ee.abs().sqrt() * ur;
            // Eq. (49)
            lambdau = su * (tu01 - neg1_pow(tpr) * tu11) / ee.abs().sqrt();
        } else if u0 >= u3 {
            // Table 1 Row 8 (line 2235)
            if firstpt && u0 != u3 {
                p[4] = -2;
                // First three u0 integrals in Eq. (47)
                phiu01 = -ellquarticreal(
                    &p,
                    -u1,
                    1.0,
                    -u2,
                    1.0,
                    -u3,
                    1.0,
                    u4,
                    -1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu0,
                    u3,
                    u0,
                );
                phiu02 = -ellquarticreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 1.0, -umi, &mut rffu0, u3, u0,
                );
                tu02 = -ellquarticreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u3, u0,
                );
                tu03 = phiu01;
                tu04 = phiu02;
                p[4] = -4;
                // Fourth u0 integral in Eq. (47)
                tu01 = -ellquarticreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu0, u3, u0,
                );
            } else if u0 == u3 {
                tu01 = 0.0;
                tu02 = 0.0;
                tu03 = 0.0;
                tu04 = 0.0;
                phiu01 = 0.0;
                phiu02 = 0.0;
            } else {
                phiu01 = tu03;
                phiu02 = tu04;
            }
            if uf != u3 {
                p[4] = -2;
                // First three uf integrals in Eq. (47)
                phiu11 = -ellquarticreal(
                    &p,
                    -u1,
                    1.0,
                    -u2,
                    1.0,
                    -u3,
                    1.0,
                    u4,
                    -1.0,
                    1.0,
                    -1.0 / uplus,
                    &mut rffu1,
                    u3,
                    uf,
                );
                phiu12 = -ellquarticreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 1.0, -umi, &mut rffu1, u3, uf,
                );
                tu12 = -ellquarticreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu1, u3, uf,
                );
                tu13 = phiu11;
                tu14 = phiu12;
                p[4] = -4;
                // Fourth uf integral in Eq. (47)
                tu11 = -ellquarticreal(
                    &p, -u1, 1.0, -u2, 1.0, -u3, 1.0, u4, -1.0, 0.0, 1.0, &mut rffu1, u3, uf,
                );
            } else {
                tu11 = 0.0;
                tu12 = 0.0;
                tu13 = 0.0;
                tu14 = 0.0;
                phiu11 = 0.0;
                phiu12 = 0.0;
            }
            // Eq. (47) for u0 piece
            tu0 = (umi - 1.0 / uplus) * tu01 + (umi * umi - 1.0 / (uplus * uplus)) * tu02
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu03
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu04;
            // Eq. (47) for uf piece
            tu1 = (umi - 1.0 / uplus) * tu11 + (umi * umi - 1.0 / (uplus * uplus)) * tu12
                - (2.0 * a * (a - l) + a * a / uplus + 1.0 / (uplus * uplus * uplus)) * tu13
                + (2.0 * a * (a - l) + a * a * umi + umi * umi * umi) * tu14;
            // Eq. (47)
            tu = su * (tu0 - neg1_pow(tpr) * tu1) / ee.abs().sqrt() * ur;
            // Eq. (48) for u0 piece
            phiu0 = -(l / uplus + 2.0 * (a - l)) * phiu01 + (l * umi + 2.0 * (a - l)) * phiu02;
            // Eq. (48) for uf piece
            phiu1 = -(l / uplus + 2.0 * (a - l)) * phiu11 + (l * umi + 2.0 * (a - l)) * phiu12;
            // Eq. (48)
            phiu = su * (phiu0 - neg1_pow(tpr) * phiu1) / ee.abs().sqrt() * ur;
            // Eq. (49) for u term
            lambdau = su * (tu01 - neg1_pow(tpr) * tu11) / ee.abs().sqrt();
        }
    }
    if ncase > 4 {
        // Find roots of biquadratic M(mu) (lines 2297-2349)
        let aa_ql2 = a * a - ql2;
        let yy = -0.5 * (aa_ql2 + fsgn1(aa_ql2) * (aa_ql2 * aa_ql2 + 4.0 * q2 * a * a).sqrt());
        let (mut mneg, mut mpos);
        if aa_ql2 < 0.0 {
            mneg = -yy / a / a;
            mpos = q2 / yy;
        } else {
            mneg = q2 / yy;
            mpos = -yy / a / a;
        }
        // Protect against rounding error in mpos
        if mpos > 1.0 {
            mpos = 1.0;
        }
        muplus = mpos.sqrt();
        // NOTE: This formula uses a slightly different prescription for
        // splitting up the mu integral. All integrals are with respect to
        // muplus, but the procedure of splitting into coefficients and
        // finding them by writing out specific cases is the same.
        a1 = sm;
        a2 = sm * neg1_pow(tpm + 1);
        a3 = 2.0 * (((2.0 * f64::from(tpm) - sm + 1.0) / 4.0).trunc());
        if mneg < 0.0 {
            // Protect against rounding errors in roots
            if muplus < mu0 {
                muplus = mu0;
                mpos = muplus * muplus;
            }
            // This is the symmetric roots case, where the orbit can cross the
            // equatorial plane. Calculate phi, t mu component integrals.
            let out = calcphitmusym(
                a, mneg, mpos, mu0, muf, muplus, rffmu1, rffmu2, rffmu3, rdc, rjc, firstpt,
            );
            rdc = out.rdc;
            rjc = out.rjc;
            if firstpt {
                tmu1 = out.tmu0;
                tmu3 = out.tmum;
                phimu1 = out.phimu0;
                phimu3 = out.phimum;
            }
            tmu2 = out.tmuf;
            phimu2 = out.phimuf;
            // Eq. (45)
            tmu = a * a * mneg * iu + (a1 * tmu1 + a2 * tmu2 + a3 * tmu3);
            // Eq. (46)
            phimu = -l * iu + l * (a1 * phimu1 + a2 * phimu2 + a3 * phimu3);
        } else {
            // This is the asymmetric roots case
            if fsgn1(mu0) == -1.0 {
                muplus = -muplus;
            }
            // Protect for rounding error when mu0 is a turning point
            if muplus.abs() < mu0.abs() {
                muplus = mu0;
                mpos = muplus * muplus;
            }
            if mneg.abs() > mu0 * mu0 {
                mneg = mu0 * mu0;
            }
            // Calculate phi, t mu component integrals
            let out = calcphitmuasym(
                a, mneg, mpos, mu0, muf, muplus, rffmu1, rffmu2, rffmu3, rdc, rjc, firstpt,
            );
            rdc = out.rdc;
            rjc = out.rjc;
            if firstpt {
                tmu1 = out.tmu0;
                tmu3 = out.tmum;
                phimu1 = out.phimu0;
                phimu3 = out.phimum;
            }
            tmu2 = out.tmuf;
            phimu2 = out.phimuf;
            // Eq. (45)
            tmu = a1 * tmu1 + a2 * tmu2 + a3 * tmu3;
            // Eq. (46)
            phimu = -l * iu + l * (a1 * phimu1 + a2 * phimu2 + a3 * phimu3);
        }
    }
    // Eq. (49)
    let lambda = lambdau + tmu;
    if l == 0.0 {
        phiu = phiu - fsgn1(a) * pi * f64::from(tpm);
    }
    PhitimeOutput {
        phimu,
        tmu,
        phiu,
        tu,
        lambda,
        state: PhitimeState {
            rdc,
            rjc,
            tu01,
            tu02,
            tu03,
            tu04,
            tmu1,
            tmu3,
            phimu1,
            phimu3,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Equatorial photon with `a = 0` and `q2 > 0`: the mu integrals are
    /// elementary (lines 1941-1965) and the u integrals vanish for
    /// `ncase = 0`. For `mu0 = muf = 0` the mu terms collapse to
    /// `phimu ~ 0` and `lambda = tmu + lambdau = 0`; everything must be
    /// finite and the state fields not written by this path must pass
    /// through unchanged on `firstpt = false`.
    #[test]
    fn equatorial_a_zero_is_finite() {
        let a = 0.0;
        let l = 1.5;
        let l2 = l * l;
        let q2 = 1.0;
        let state = PhitimeState::default();
        let out = geophitime(
            0.3, 0.4, 0.0, 0.0, a, l, l2, q2, 0, 0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0, &state, true,
        );
        assert!(out.lambda.is_finite());
        assert!(out.tu.is_finite());
        assert!(out.tmu.is_finite());
        assert!(out.phiu.is_finite());
        assert!(out.phimu.is_finite());
        // ncase = 0: all u terms are exactly zero, and the a = 0 mu terms
        // cancel for mu0 = muf = 0.
        assert_eq!(out.tu, 0.0);
        assert_eq!(out.tmu, 0.0);
        assert_eq!(out.lambda, 0.0);
        assert!(out.phimu.abs() < 1e-12);
        assert_eq!(out.phiu, 0.0);

        // On a later point (firstpt = false) the cached state fields that this
        // path never writes come from the caller.
        let cached = PhitimeState {
            rdc: 1.25,
            rjc: -3.5,
            tu01: 0.1,
            tu02: 0.2,
            tu03: 0.3,
            tu04: 0.4,
            tmu1: 0.5,
            tmu3: 0.6,
            phimu1: 0.7,
            phimu3: 0.8,
        };
        let out2 = geophitime(
            0.3, 0.4, 0.0, 0.0, a, l, l2, q2, 0, 0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0, &cached, false,
        );
        assert_eq!(out2.state.rdc, cached.rdc);
        assert_eq!(out2.state.rjc, cached.rjc);
        assert_eq!(out2.state.tu01, cached.tu01);
        assert_eq!(out2.state.tu02, cached.tu02);
        assert_eq!(out2.state.tu04, cached.tu04);
        assert_eq!(out2.state.tmu1, cached.tmu1);
        assert_eq!(out2.state.tmu3, cached.tmu3);
        // phimu1 and phimu3 are recomputed by the a = 0 branch.
        assert!(out2.state.phimu1.is_finite());
        assert!((out2.state.phimu3 - 2.0 * (-1.0_f64).acos()).abs() < 1e-12);
    }

    /// The `ncase = 4` path (q2 = 0, l = a) is elementary, so it can be
    /// checked before the elliptic-integral helpers of the radial motion are
    /// ported. As `uf -> u0` both the time and azimuth increments tend to
    /// zero.
    #[test]
    fn small_step_time_and_azimuth_vanish_ncase4() {
        let a = 0.5;
        let l = 0.5;
        let q2 = 0.0;
        let u0 = 0.2;
        let uf = u0 * (1.0 + 1.0e-9);
        let state = PhitimeState::default();
        let out = geophitime(
            u0,
            uf,
            0.3,
            0.3,
            a,
            l,
            l * l,
            q2,
            0,
            0,
            -1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            4,
            &state,
            true,
        );
        assert!(out.tu.abs() < 1.0e-6, "tu = {}", out.tu);
        assert!(out.phiu.abs() < 1.0e-6, "phiu = {}", out.phiu);
        assert!(out.lambda.abs() < 1.0e-6, "lambda = {}", out.lambda);
        assert!(out.tu.is_finite() && out.phiu.is_finite() && out.lambda.is_finite());
    }

    /// Small-step consistency in a quartic (`ncase = 7`) case. This exercises
    /// `ellquarticreal`, which is still `todo!()` in the elliptic module
    /// (implemented by the radial port), so the test is ignored until that
    /// lands.
    #[test]
    #[ignore = "needs elliptic::ellquarticreal (still todo!() in the elliptic port)"]
    fn small_step_time_and_azimuth_vanish_ncase7() {
        let a = 0.9;
        let l = 2.0;
        let l2 = l * l;
        let q2 = 4.0;
        let u0 = 0.1;
        let uf = u0 * (1.0 + 1.0e-9);
        // Root data would normally come from GEOMU; placeholders keep this
        // test compiling until the full call chain is available.
        let state = PhitimeState::default();
        let out = geophitime(
            u0, uf, 0.6, 0.6, a, l, l2, q2, 0, 0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 7, &state, true,
        );
        assert!(out.tu.is_finite() && out.phiu.is_finite() && out.lambda.is_finite());
    }
}
