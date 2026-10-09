//! Special functions used by geokerr: Gauss-Legendre quadrature, polynomial
//! root finding (ZROOTS/LAGUER) and Jacobi elliptic functions (SNCNDN).
//!
//! Translation target: `geokerr_wrapper.f` routines GAULEG, ZROOTS, LAGUER,
//! SNCNDN.

use num_complex::Complex64;

/// Gauss-Legendre abscissas and weights on [x1, x2] with n points
/// (upstream `gauleg`). Returns `(x, w)`.
///
/// Subroutine to calculate the abscissas and weights for the Gauss-Legendre
/// quadrature (Numerical Recipes; algorithm of G. B. Rybicki). The returned
/// `x`/`w` vectors are ordered as the upstream arrays `X(1..n)`, `W(1..n)`.
pub fn gauleg(x1: f64, x2: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    const PI: f64 = std::f64::consts::PI;
    const EPS: f64 = 3e-14;

    let m = n.div_ceil(2);
    let xm = 0.5 * (x2 + x1);
    let xl = 0.5 * (x2 - x1);
    let mut x = vec![0.0f64; n];
    let mut w = vec![0.0f64; n];
    for i in 1..=m {
        let mut z = (PI * (i as f64 - 0.25) / (n as f64 + 0.5)).cos();
        let pp;
        // Newton iteration on the Legendre polynomial (upstream GOTO 1).
        loop {
            let mut p1 = 1.0f64;
            let mut p2 = 0.0f64;
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j as f64 - 1.0) * z * p2 - (j as f64 - 1.0) * p3) / j as f64;
            }
            let pp_local = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / pp_local;
            if (z - z1).abs() <= EPS {
                pp = pp_local;
                break;
            }
        }
        x[i - 1] = xm - xl * z;
        x[n - i] = xm + xl * z;
        w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - i] = w[i - 1];
    }
    (x, w)
}

/// Find all roots of the complex polynomial `a[0..=m]` (upstream `zroots`).
/// Returns the `m` roots in `roots[0..m]`.
///
/// Given the degree `m` and the `m+1` complex coefficients `a` of the
/// polynomial (with `a[0]` the constant term), returns all `m` roots in the
/// complex vector `roots`. `polish` requests polishing of the roots by
/// Laguerre's method with the undeflated coefficients (Press et al. 1992).
pub fn zroots(a: &[Complex64], m: usize, polish: bool) -> Vec<Complex64> {
    const EPS: f64 = 1e-6;
    const MAXM: usize = 7;

    if m > MAXM - 1 {
        // Upstream: WRITE(6,*) 'M too large in ZROOTS' and continues.
        eprintln!("M too large in ZROOTS");
    }
    // Copy of coefficients for successive deflation.
    let mut ad: Vec<Complex64> = a[..=m].to_vec();
    // Upstream also sets XSUM=0 here; it is never used (the accumulation
    // lines are commented out).
    let mut xstart = Complex64::new(0.0, 0.0);
    // If ncc=1, the previous root is a complex conjugate of another root.
    let mut ncc: i32 = 0;
    let mut roots = vec![Complex64::new(0.0, 0.0); m];
    // Loop over each root to be found.
    for j in (1..=m).rev() {
        // Start at zero to favour convergence to the smallest remaining
        // root, or if the previous root was complex, start from the
        // upstream DCMPLX(RE(ROOT),-RE(ROOT)) (kept verbatim, including its
        // use of the real part rather than the imaginary part).
        let mut x = xstart;
        if j < m {
            if roots[j].im != 0.0 && ncc == 0 {
                xstart = Complex64::new(roots[j].re, -roots[j].re);
                // Since we have chosen the second root to start at its
                // complex conjugate, do not use the conjugate again.
                ncc = 1;
            } else {
                xstart = Complex64::new(0.0, 0.0);
                ncc = 0;
            }
        } else {
            xstart = Complex64::new(0.0, 0.0);
            ncc = 0;
        }
        // Find the root.
        x = laguer(&ad, j, x, 0);
        if x.im.abs() <= 2.0 * EPS * EPS * x.re.abs() {
            x = Complex64::new(x.re, 0.0);
        }
        roots[j - 1] = x;
        // Forward deflation.
        let mut b = ad[j];
        for jj in (1..=j).rev() {
            let c = ad[jj - 1];
            ad[jj - 1] = b;
            b = x * b + c;
        }
    }
    if polish {
        // Polish the roots using the undeflated coefficients.
        for j in 1..=m {
            roots[j - 1] = laguer(a, m, roots[j - 1], 0);
        }
    }
    // Sort roots by their real parts by straight insertion.
    for j in 2..=m {
        let x = roots[j - 1];
        let mut i = 0usize;
        for ii in (1..=j - 1).rev() {
            if roots[ii - 1].re <= x.re {
                i = ii;
                break;
            }
            roots[ii] = roots[ii - 1];
        }
        roots[i] = x;
    }
    roots
}

/// Laguerre's method iteration for one root (upstream `laguer`).
///
/// Finds one root of the polynomial `a[0..=m]` starting from `x`.
/// (`its`, the Fortran ITS iteration counter output, is accepted for
/// signature compatibility and unused.)
pub fn laguer(a: &[Complex64], m: usize, x: Complex64, its: u32) -> Complex64 {
    let _ = its;
    const MR: usize = 8;
    const MT: usize = 10;
    const MAXIT: usize = MT * MR; // 100
    const EPSS: f64 = 1e-15;
    const FRAC: [f64; MR] = [0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0];

    // Loop over iterations up to the allowed maximum.
    let mut x = x;
    for iter in 1..=MAXIT {
        let mut b = a[m];
        let mut err = b.norm();
        let mut d = Complex64::new(0.0, 0.0);
        let mut f = Complex64::new(0.0, 0.0);
        let abx = x.norm();
        for j in (1..=m).rev() {
            // Efficient computation of the polynomial and its first two
            // derivatives.
            f = x * f + d;
            d = x * d + b;
            b = x * b + a[j - 1];
            err = b.norm() + abx * err;
        }
        err *= EPSS;

        if b.norm() <= err {
            // Special case: we are on the root.
            return x;
        } else {
            // The generic case; use Laguerre's formula.
            let g = d / b;
            let g2 = g * g;
            let h = g2 - 2.0 * f / b;
            let sq = ((m as f64 - 1.0) * (m as f64 * h - g2)).sqrt();
            let mut gp = g + sq;
            let gm = g - sq;
            let abp = gp.norm();
            let abm = gm.norm();
            if abp < abm {
                gp = gm;
            }
            let dx = if abp.max(abm) > 0.0 {
                m as f64 / gp
            } else {
                Complex64::new((1.0 + abx).ln(), iter as f64).exp()
            };
            let x1 = x - dx;
            // Check if we have converged.
            if x == x1 {
                return x;
            }
            if iter % MT != 0 {
                x = x1;
            } else {
                x -= dx * FRAC[iter / MT - 1];
            }
        }
    }
    // Upstream: WRITE(6,*) 'Too many iterations'
    eprintln!("Too many iterations");
    x
}

/// Jacobi elliptic functions sn(u, k), cn(u, k), dn(u, k) with
/// `emmc = 1 - k^2` (upstream `sncndn`). Returns `(sn, cn, dn)`.
///
/// Given the argument u and `emmc = 1 - k^2`, computes sn(u,k), cn(u,k),
/// dn(u,k) (Press et al. 1992).
pub fn sncndn(uu: f64, emmc: f64) -> (f64, f64, f64) {
    const CA: f64 = 3e-8;

    let mut em = [0.0f64; 13];
    let mut en = [0.0f64; 13];
    let mut emc = emmc;
    let mut u = uu;
    let sn;
    let cn;
    let mut dn;
    if emc != 0.0 {
        let bo = emc < 0.0;
        let mut d = 0.0f64;
        if bo {
            d = 1.0 - emc;
            emc = -emc / d;
            d = d.sqrt();
            u *= d;
        }
        let mut a = 1.0f64;
        dn = 1.0f64;
        let mut c = 0.0f64;
        let mut l = 13usize;
        for i in 1..=13 {
            l = i;
            em[i - 1] = a;
            emc = emc.sqrt();
            en[i - 1] = emc;
            c = 0.5 * (a + emc);
            if (a - emc).abs() <= CA * a {
                break;
            }
            emc *= a;
            a = c;
        }
        u *= c;
        let mut snv = u.sin();
        let mut cnv = u.cos();
        if snv != 0.0 {
            a = cnv / snv;
            c *= a;
            for ii in (1..=l).rev() {
                let b = em[ii - 1];
                a *= c;
                c *= dn;
                dn = (en[ii - 1] + a) / (b + a);
                a = c / b;
            }
            a = 1.0 / (c * c + 1.0).sqrt();
            if snv < 0.0 {
                snv = -a;
            } else {
                snv = a;
            }
            cnv = c * snv;
        }
        if bo {
            std::mem::swap(&mut dn, &mut cnv);
            snv /= d;
        }
        sn = snv;
        cn = cnv;
    } else {
        cn = 1.0 / u.cosh();
        dn = cn;
        sn = u.tanh();
    }
    (sn, cn, dn)
}

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
    fn sncndn_at_zero_and_hyperbolic_limit() {
        // u = 0, k = sqrt(1/2) (emmc = 0.5): sn=0, cn=dn=1.
        let (sn, cn, dn) = sncndn(0.0, 0.5);
        assert_eq!(sn, 0.0);
        assert_eq!(cn, 1.0);
        assert_eq!(dn, 1.0);

        // emmc = 0 <=> k = 1: sn=tanh, cn=dn=sech.
        for &u in &[-2.0f64, -0.5, 0.25, 1.0, 2.0] {
            let (sn, cn, dn) = sncndn(u, 0.0);
            assert_close(sn, u.tanh(), 1e-13, "sn(u,1)");
            assert_close(cn, 1.0 / u.cosh(), 1e-13, "cn(u,1)");
            assert_close(dn, 1.0 / u.cosh(), 1e-13, "dn(u,1)");
        }
    }

    #[test]
    fn sncndn_circular_limit() {
        // emmc = 1 <=> k = 0: sn=sin, cn=cos, dn=1.
        for &u in &[-2.0f64, -0.5, 0.25, 1.0, 2.0] {
            let (sn, cn, dn) = sncndn(u, 1.0);
            assert_close(sn, u.sin(), 1e-13, "sn(u,0)");
            assert_close(cn, u.cos(), 1e-13, "cn(u,0)");
            assert_close(dn, 1.0, 1e-13, "dn(u,0)");
        }
    }

    #[test]
    fn gauleg_weights_sum_and_symmetry() {
        let (x, w) = gauleg(-1.0, 1.0, 8);
        assert_eq!(x.len(), 8);
        assert_eq!(w.len(), 8);
        let sum: f64 = w.iter().sum();
        assert!((sum - 2.0).abs() < 1e-13, "sum(w) = {sum}");
        for i in 0..8 {
            assert!((x[i] + x[7 - i]).abs() < 1e-13, "x symmetry at {i}");
            assert!((w[i] - w[7 - i]).abs() < 1e-13, "w symmetry at {i}");
        }
        // n=8 integrates polynomials up to degree 15: \int x^2 = 2/3.
        let m2: f64 = x.iter().zip(&w).map(|(xi, wi)| wi * xi * xi).sum();
        assert!((m2 - 2.0 / 3.0).abs() < 1e-13, "int x^2 = {m2}");
    }

    #[test]
    fn zroots_finds_quartic_roots() {
        // (z-1)(z-2)(z-3)(z-4) = z^4 - 10 z^3 + 35 z^2 - 50 z + 24.
        let a = [
            Complex64::new(24.0, 0.0),
            Complex64::new(-50.0, 0.0),
            Complex64::new(35.0, 0.0),
            Complex64::new(-10.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        let roots = zroots(&a, 4, true);
        assert_eq!(roots.len(), 4);
        let mut re: Vec<f64> = roots.iter().map(|r| r.re).collect();
        for r in &roots {
            assert!(r.im.abs() < 1e-10, "imaginary part of root {r}");
        }
        re.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, r) in re.iter().enumerate() {
            assert!(
                (r - (i as f64 + 1.0)).abs() < 1e-10,
                "root {i}: got {r}, want {}",
                i + 1
            );
        }
    }

    #[test]
    fn zroots_finds_complex_pair() {
        // (z-1)((z-2)^2+9) = z^3 - 5 z^2 + 17 z - 13.
        let a = [
            Complex64::new(-13.0, 0.0),
            Complex64::new(17.0, 0.0),
            Complex64::new(-5.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        let roots = zroots(&a, 3, true);
        let want = [
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 3.0),
            Complex64::new(2.0, -3.0),
        ];
        // Roots come back sorted by real part; match greedily.
        for w in want {
            let best = roots
                .iter()
                .map(|r| (r - w).norm())
                .fold(f64::INFINITY, f64::min);
            assert!(best < 1e-9, "no root near {w} (best {best})");
        }
    }
}
