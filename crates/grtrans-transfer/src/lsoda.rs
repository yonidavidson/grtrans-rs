//! Adaptive ODE integration for the transfer equation (upstream `lsoda`).
//!
//! ## Documented substitution
//!
//! Upstream integrates the polarized transfer equation with ODEPACK's LSODA
//! (rtol=1e-6, atol=1e-8, hmax=0.1 in affine units) on the RHS
//! `dI/ds = j(s) - K(s) I` with piecewise-linear coefficient interpolation.
//! Porting the 28k-line F77 ODEPACK sources is out of scope for this port,
//! so this module provides a Dormand-Prince 5(4) adaptive solver with the
//! same tolerances and step limit. The solver is validated against the
//! upstream regression problems (FFJET etc., see tests/reference) and
//! cross-checked against the exact `delo`/`formal` schemes; the deviation
//! is documented in docs/PORTING_MATRIX.md and docs/VALIDATION_PLAN.md.

use crate::opacity_matrix;

/// RHS: `dI/ds = j(s) - K(s) I` with linear interpolation of `j` and `K`
/// between tabulated samples (upstream `radtrans_aux` + `radtrans_rhs_form`).
fn rhs(ss: &[f64], jj: &[[f64; 4]], kk: &[[[f64; 4]; 4]], lam: f64, y: &[f64; 4]) -> [f64; 4] {
    let n = ss.len();
    // locate interval (ss ascending)
    let mut k = 0usize;
    if lam <= ss[0] {
        k = 0;
    } else if lam >= ss[n - 1] {
        k = n - 2;
    } else {
        // binary search for ss[k] <= lam < ss[k+1]
        let mut lo = 0usize;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if lam >= ss[mid] {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        k = lo;
    }
    let w = (lam - ss[k]) / (ss[k + 1] - ss[k]);
    let mut jv = [0.0f64; 4];
    let mut km = [[0.0f64; 4]; 4];
    for i in 0..4 {
        jv[i] = (1.0 - w) * jj[k][i] + w * jj[k + 1][i];
        for q in 0..4 {
            km[i][q] = (1.0 - w) * kk[k][i][q] + w * kk[k + 1][i][q];
        }
    }
    let mut d = [0.0f64; 4];
    for i in 0..4 {
        let mut acc = 0.0;
        for q in 0..4 {
            acc += km[i][q] * y[q];
        }
        d[i] = jv[i] - acc;
    }
    d
}

/// Dormand-Prince 5(4) adaptive integration of the transfer equation.
///
/// `s` is the affine-parameter array as passed to the other integrators
/// (decreasing along the ray); the integration internally reverses it to be
/// increasing, exactly as upstream does before calling LSODA. `tau` is the
/// positive optical depth array; the integration range is trimmed to the
/// region with nonzero emission and tau <= MAX_TAU (10), as upstream.
///
/// Returns `(intensity, nptsout)` with the intensity (4 x npts) at the
/// original sample points and the upstream output index.
pub fn integrate_lsoda(
    s: &[f64],
    j: &[[f64; 4]],
    a: &[[f64; 4]],
    rho: &[[f64; 3]],
    tau: &[f64],
    atol: f64,
    rtol: f64,
    hmax: f64,
) -> (Vec<[f64; 4]>, usize) {
    const MAX_TAU: f64 = 10.0;
    let npts = s.len();
    let mut intensity = vec![[0.0f64; 4]; npts];
    if npts == 1 {
        intensity[0] = j[0];
        return (intensity, 1);
    }
    // upstream: lamdex = npts if maxval(tau) <= MAX_TAU, else locate+1
    let mut lamdex = npts;
    if tau.iter().cloned().fold(f64::NEG_INFINITY, f64::max) > MAX_TAU {
        // locate(tau, MAX_TAU): 1-based bracket index on an ascending array
        let mut jl = 0isize;
        let mut ju = npts as isize + 1;
        while ju - jl > 1 {
            let jm = (ju + jl) / 2;
            if MAX_TAU >= tau[(jm - 1) as usize] {
                jl = jm;
            } else {
                ju = jm;
            }
        }
        lamdex = (jl + 1) as usize;
    }
    // only use the parts of the ray where emissivity is non-zero
    let mut i1 = 0usize; // 0-based
    if j[0][0] == 0.0 {
        for ii in 1..npts {
            if j[ii][0] != 0.0 {
                i1 = ii;
                break;
            }
        }
    }
    let mut i2 = lamdex - 1; // 0-based inclusive
    if j[lamdex - 1][0] == 0.0 {
        for ii in 1..lamdex {
            if j[lamdex - 1 - ii][0] != 0.0 {
                i2 = lamdex - 1 - ii;
                break;
            }
        }
    }
    let nptsout = i2 + 1;
    // integration range (0-based, inclusive): i1..=i2; upstream s(i1:i2)
    let m = i2 - i1 + 1;
    // upstream reverses: ss(1:npts)=s0(npts:1:-1); s(i1:i2)=s(i2:i1:-1)
    let ss: Vec<f64> = (0..m).map(|k| s[i2 - k]).collect();
    let jj: Vec<[f64; 4]> = (0..m).map(|k| j[i2 - k]).collect();
    let karr: Vec<[[f64; 4]; 4]> = (0..m)
        .map(|k| opacity_matrix(&a[i2 - k], &rho[i2 - k]))
        .collect();
    if m == 1 {
        intensity[i2] = jj[0];
        return (intensity, nptsout);
    }
    // output at the reversed points, then map back
    let mut yout = vec![[0.0f64; 4]; m];
    let mut y = [0.0f64; 4];
    yout[0] = y;
    let mut x = ss[0];
    let mut h = (ss[1] - ss[0]).min(hmax);
    let mut out_idx = 1usize;
    let mut steps = 0usize;
    let max_steps = 10_000_000usize;
    // Dormand-Prince 5(4) coefficients
    let c2 = 1.0 / 5.0;
    let c3 = 3.0 / 10.0;
    let c4 = 4.0 / 5.0;
    let c5 = 8.0 / 9.0;
    while out_idx < m && steps < max_steps {
        steps += 1;
        let x_target = ss[out_idx];
        if x + h > x_target {
            h = x_target - x;
        }
        if h < 1e-300 {
            h = (x_target - x).max(1e-300);
        }
        let k1 = rhs(&ss, &jj, &karr, x, &y);
        let mut yt = [0.0f64; 4];
        for i in 0..4 {
            yt[i] = y[i] + h * (c2 * k1[i]);
        }
        let k2 = rhs(&ss, &jj, &karr, x + c2 * h, &yt);
        for i in 0..4 {
            yt[i] = y[i] + h * (3.0 / 40.0 * k1[i] + 9.0 / 40.0 * k2[i]);
        }
        let k3 = rhs(&ss, &jj, &karr, x + c3 * h, &yt);
        for i in 0..4 {
            yt[i] = y[i] + h * (44.0 / 45.0 * k1[i] - 56.0 / 15.0 * k2[i] + 32.0 / 9.0 * k3[i]);
        }
        let k4 = rhs(&ss, &jj, &karr, x + c4 * h, &yt);
        for i in 0..4 {
            yt[i] = y[i]
                + h * (19372.0 / 6561.0 * k1[i] - 25360.0 / 2187.0 * k2[i]
                    + 64448.0 / 6561.0 * k3[i]
                    - 212.0 / 729.0 * k4[i]);
        }
        let k5 = rhs(&ss, &jj, &karr, x + c5 * h, &yt);
        for i in 0..4 {
            yt[i] = y[i]
                + h * (9017.0 / 3168.0 * k1[i] - 355.0 / 33.0 * k2[i]
                    + 46732.0 / 5247.0 * k3[i]
                    + 49.0 / 176.0 * k4[i]
                    - 5103.0 / 18656.0 * k5[i]);
        }
        let k6 = rhs(&ss, &jj, &karr, x + h, &yt);
        let mut y5 = [0.0f64; 4];
        for i in 0..4 {
            y5[i] = y[i]
                + h * (35.0 / 384.0 * k1[i] + 500.0 / 1113.0 * k3[i] + 125.0 / 192.0 * k4[i]
                    - 2187.0 / 6784.0 * k5[i]
                    + 11.0 / 84.0 * k6[i]);
        }
        let k7 = rhs(&ss, &jj, &karr, x + h, &y5);
        // error estimate (difference of 5th and 4th order solutions)
        let mut err = 0.0f64;
        for i in 0..4 {
            let e = h
                * (71.0 / 57600.0 * k1[i] - 71.0 / 16695.0 * k3[i] + 71.0 / 1920.0 * k4[i]
                    - 17253.0 / 339200.0 * k5[i]
                    + 22.0 / 525.0 * k6[i]
                    - 1.0 / 40.0 * k7[i]);
            let sc = atol + rtol * y[i].abs().max(y5[i].abs());
            err = err.max((e / sc).abs());
        }
        if err <= 1.0 {
            // accept
            x += h;
            y = y5;
            if (x - ss[out_idx]).abs() < 1e-15 * ss[out_idx].abs().max(1.0) || x >= ss[out_idx] {
                yout[out_idx] = y;
                out_idx += 1;
            }
        }
        // adapt step
        let fac = if err == 0.0 {
            5.0
        } else {
            (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
        };
        h *= fac;
        h = h.min(hmax);
        if h <= 0.0 {
            h = 1e-12;
        }
    }
    if out_idx < m {
        // fill any remaining outputs by holding the last value (should not
        // happen with max_steps generous)
        for k in out_idx..m {
            yout[k] = y;
        }
    }
    // Map back following the upstream array-section semantics: the
    // solution at s(i1 + k) (section order, i.e. reversed) is yout[k], so
    // the observer value (largest s) lands at index i2, as upstream.
    for k in 0..m {
        intensity[i1 + k] = yout[k];
    }
    (intensity, nptsout)
}
