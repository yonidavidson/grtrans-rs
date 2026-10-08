//! Mathematical utilities: cumulative trapezoidal integration, dot products
//! and Brent root finding.
//!
//! Direct translation of `math.f90` (upstream GRTRANS).

/// Cumulative trapezoidal integral with the same convention as upstream
/// `tsum`:
///
/// ```text
/// tsum[i] = sum_{j=1..i-1} (x[j+1]-x[j]) * (y[j+1]+y[j])/2
/// ```
///
/// i.e. `tsum[0] = 0` and `tsum[k]` is the integral from `x[0]` to `x[k]`
/// computed by the trapezoid rule. The return value has the same length as
/// the inputs.
pub fn tsum(x: &[f64], y: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), y.len());
    let n = x.len();
    let mut out = vec![0.0f64; n];
    // upstream: tsumint = [0, dx_i * yavg_i], then cumulative sum
    let mut acc = 0.0f64;
    for i in 1..n {
        let xdif = x[i] - x[i - 1];
        let yavg = (y[i - 1] + y[i]) / 2.0;
        acc += xdif * yavg;
        out[i] = acc;
    }
    out
}

/// Cumulative sum, as upstream `cum_sum` (each output element is the sum of
/// all preceding inputs up to and including the current one).
pub fn cum_sum(x: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0f64; x.len()];
    let mut acc = 0.0f64;
    for (i, v) in x.iter().enumerate() {
        acc += *v;
        out[i] = acc;
    }
    out
}

/// Row-wise dot product of two row-major `nrows x ncols` arrays, matching
/// upstream `dot_product_arr` (result length `nrows`).
pub fn dot_product_arr(a: &[f64], b: &[f64], nrows: usize, ncols: usize) -> Vec<f64> {
    assert_eq!(a.len(), nrows * ncols);
    assert_eq!(b.len(), nrows * ncols);
    let mut dot = vec![0.0f64; nrows];
    for i in 0..nrows {
        let mut acc = 0.0f64;
        for j in 0..ncols {
            acc += a[i * ncols + j] * b[i * ncols + j];
        }
        dot[i] = acc;
    }
    dot
}

/// Brent root finder, direct translation of upstream `zbrent`
/// (Numerical-Recipes-style algorithm).
///
/// Panics if the root is not bracketed, mirroring upstream's `stop`.
pub fn zbrent<F: Fn(f64) -> f64>(func: F, x1: f64, x2: f64, tol: f64) -> f64 {
    const ITMAX: usize = 100;
    let eps = f64::EPSILON;
    let mut a = x1;
    let mut b = x2;
    let mut fa = func(a);
    let mut fb = func(b);
    if (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0) {
        panic!("zbrent: root must be bracketed: x1={x1}, x2={x2}, fa={fa}, fb={fb}");
    }
    let mut c = b;
    let mut fc = fb;
    #[allow(unused_assignments)]
    // Fortran leaves d,e undefined on the first iteration, but the branch
    // that assigns them (fc=fb makes the sign test true) always executes
    // before any read, except in the fb==0 case where the function returns
    // first. Initializing to zero is therefore behaviour-preserving.
    let mut d = 0.0f64; // assigned before any read; see comment
    let mut e = 0.0f64;
    for _ in 0..ITMAX {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        } else {
            d = e;
        }
        if fc.abs() < fb.abs() {
            // Fortran: a=b; b=c; c=a; fa=fb; fb=fc; fc=fa (sequential)
            let old_b = b;
            let old_c = c;
            let old_fb = fb;
            let old_fc = fc;
            a = old_b;
            b = old_c;
            c = a;
            fa = old_fb;
            fb = old_fc;
            fc = fa;
        }
        let tol1 = 2.0 * eps * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return b;
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (p, q);
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            let mut p = p;
            let mut qq = q;
            if p > 0.0 {
                qq = -qq;
            }
            p = p.abs();
            if 2.0 * p < (3.0 * xm * qq - (tol1 * qq).abs()).min((e * qq).abs()) {
                e = d;
                d = p / qq;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol1 { d } else { tol1.copysign(xm) };
        fb = func(b);
    }
    // upstream prints a warning and returns the last b
    b
}

/// Vectorized Brent root finder, direct translation of upstream
/// `zbrent_array`. Each element is solved independently; the Fortran
/// "converged elements freeze" semantics are reproduced exactly.
pub fn zbrent_array<F: Fn(&[f64]) -> Vec<f64>>(
    func: F,
    x1: &[f64],
    x2: &[f64],
    tol: f64,
) -> Vec<f64> {
    const ITMAX: usize = 100;
    let eps = f64::EPSILON;
    let n = x1.len();
    assert_eq!(x2.len(), n);
    let mut a = x1.to_vec();
    let mut b = x2.to_vec();
    let mut fa = func(&a);
    let mut fb = func(&b);
    let mut c = b.clone();
    let mut fc = fb.clone();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    let mut zbrent = b.clone();
    let mut converge = vec![false; n];
    for _ in 0..ITMAX {
        for i in 0..n {
            if converge[i] {
                continue;
            }
            if (fb[i] > 0.0 && fc[i] > 0.0) || (fb[i] < 0.0 && fc[i] < 0.0) {
                c[i] = a[i];
                fc[i] = fa[i];
                d[i] = b[i] - a[i];
                e[i] = d[i];
            }
            if fc[i].abs() < fb[i].abs() {
                let old_b = b[i];
                let old_c = c[i];
                let old_fb = fb[i];
                let old_fc = fc[i];
                a[i] = old_b;
                b[i] = old_c;
                c[i] = a[i];
                fa[i] = old_fb;
                fb[i] = old_fc;
                fc[i] = fa[i];
            }
        }
        let mut all = true;
        for i in 0..n {
            if converge[i] {
                continue;
            }
            let tol1 = 2.0 * eps * b[i].abs() + 0.5 * tol;
            let xm = 0.5 * (c[i] - b[i]);
            converge[i] = xm.abs() <= tol1 || fb[i] == 0.0;
            if converge[i] {
                zbrent[i] = b[i];
            } else {
                all = false;
            }
        }
        if all {
            return zbrent;
        }
        // update non-converged elements
        let mut new_b_eval: Vec<usize> = Vec::new();
        for i in 0..n {
            if converge[i] {
                continue;
            }
            let xm = 0.5 * (c[i] - b[i]);
            let tol1 = 2.0 * eps * b[i].abs() + 0.5 * tol;
            if e[i].abs() >= tol1 && fa[i].abs() > fb[i].abs() {
                let s = fb[i] / fa[i];
                let (mut p, mut q);
                if a[i] == c[i] {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    let qq = fa[i] / fc[i];
                    let r = fb[i] / fc[i];
                    p = s * (2.0 * xm * qq * (qq - r) - (b[i] - a[i]) * (r - 1.0));
                    q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if p > 0.0 {
                    q = -q;
                }
                p = p.abs();
                if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e[i] * q).abs()) {
                    e[i] = d[i];
                    d[i] = p / q;
                } else {
                    d[i] = xm;
                    e[i] = d[i];
                }
            } else {
                d[i] = xm;
                e[i] = d[i];
            }
            a[i] = b[i];
            fa[i] = fb[i];
            b[i] += if d[i].abs() > tol1 {
                d[i]
            } else {
                tol1.copysign(xm)
            };
            new_b_eval.push(i);
        }
        let bvals: Vec<f64> = new_b_eval.iter().map(|&i| b[i]).collect();
        let fvals = func(&bvals);
        for (k, &i) in new_b_eval.iter().enumerate() {
            fb[i] = fvals[k];
        }
    }
    zbrent
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsum_matches_trapezoid() {
        let x = [0.0, 1.0, 3.0];
        let y = [0.0, 2.0, 2.0];
        let s = tsum(&x, &y);
        assert_eq!(s[0], 0.0);
        assert_eq!(s[1], 1.0); // trapezoid over [0,1]: 1.0
        assert_eq!(s[2], 5.0); // + trapezoid over [1,3]: 4.0
    }

    #[test]
    fn zbrent_finds_root() {
        let f = |x: f64| x * x - 2.0;
        let r = zbrent(f, 0.0, 2.0, 1e-14);
        assert!((r - std::f64::consts::SQRT_2).abs() < 1e-13);
    }

    #[test]
    fn zbrent_array_finds_roots() {
        let f = |x: &[f64]| x.iter().map(|v| v * v - 2.0).collect();
        let r = zbrent_array(f, &[0.0, 1.0, 1.3], &[2.0, 1.5, 2.0], 1e-14);
        for v in r {
            assert!((v - std::f64::consts::SQRT_2).abs() < 1e-13);
        }
    }
}
