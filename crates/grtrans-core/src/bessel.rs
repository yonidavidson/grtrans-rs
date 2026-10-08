//! Modified Bessel functions I0, I1, K0, K1, Kn.
//!
//! Direct translation of `bessel.f90` (upstream GRTRANS; Numerical
//! Recipes-style polynomial approximations). Evaluation order in the
//! Horner scheme follows upstream exactly.

/// Polynomial evaluation with coefficients in ascending order using
/// Horner's scheme from the highest coefficient (upstream `poly`).
#[inline]
pub fn poly(x: f64, coeffs: &[f64]) -> f64 {
    if coeffs.is_empty() {
        return 0.0;
    }
    let mut p = coeffs[coeffs.len() - 1];
    for i in (0..coeffs.len() - 1).rev() {
        p = x * p + coeffs[i];
    }
    p
}

/// Modified Bessel function I0 (upstream `besseli0`).
pub fn besseli0(x: f64) -> f64 {
    const P: [f64; 7] = [
        1.00,
        3.51562290,
        3.08994240,
        1.20674920,
        0.26597320,
        0.360768e-10,
        0.45813e-20,
    ];
    const Q: [f64; 9] = [
        0.398942280,
        0.1328592e-10,
        0.225319e-20,
        -0.157565e-20,
        0.916281e-20,
        -0.2057706e-10,
        0.2635537e-10,
        -0.1647633e-10,
        0.392377e-20,
    ];
    let ax = x.abs();
    if ax < 3.75 {
        poly((x / 3.750).powi(2), &P)
    } else {
        let y = 3.750 / ax;
        (ax.exp() / ax.sqrt()) * poly(y, &Q)
    }
}

/// Modified Bessel function I1 (upstream `besseli1`).
pub fn besseli1(x: f64) -> f64 {
    const P: [f64; 7] = [
        0.50,
        0.878905940,
        0.514988690,
        0.150849340,
        0.2658733e-10,
        0.301532e-20,
        0.32411e-30,
    ];
    const Q: [f64; 9] = [
        0.398942280,
        -0.3988024e-10,
        -0.362018e-20,
        0.163801e-20,
        -0.1031555e-10,
        0.2282967e-10,
        -0.2895312e-10,
        0.1787654e-10,
        -0.420059e-20,
    ];
    let ax = x.abs();
    let val = if ax < 3.75 {
        ax * poly((x / 3.750).powi(2), &P)
    } else {
        let y = 3.750 / ax;
        (ax.exp() / ax.sqrt()) * poly(y, &Q)
    };
    if x < 0.0 {
        -val
    } else {
        val
    }
}

/// Modified Bessel function K0 (upstream `besselk0`).
pub fn besselk0(x: f64) -> f64 {
    const P: [f64; 7] = [
        -0.57721566,
        0.42278420,
        0.23069756,
        0.3488590e-1,
        0.262698e-2,
        0.10750e-3,
        0.74e-5,
    ];
    const Q: [f64; 7] = [
        1.25331414,
        -0.7832358e-1,
        0.2189568e-1,
        -0.1062446e-1,
        0.587872e-2,
        -0.251540e-2,
        0.53208e-3,
    ];
    if x <= 2.0 {
        let y = x * x / 4.0;
        (-(x / 2.0).ln()) * besseli0(x) + poly(y, &P)
    } else {
        let y = 2.0 / x;
        (x.exp().recip() / x.sqrt()) * poly(y, &Q)
    }
}

/// Modified Bessel function K1 (upstream `besselk1`).
pub fn besselk1(x: f64) -> f64 {
    const P: [f64; 7] = [
        1.0,
        0.15443144,
        -0.67278579,
        -0.18156897,
        -0.1919402e-1,
        -0.110404e-2,
        -0.4686e-4,
    ];
    const Q: [f64; 7] = [
        1.25331414,
        0.23498619,
        -0.3655620e-1,
        0.1504268e-1,
        -0.780353e-2,
        0.325614e-2,
        -0.68245e-3,
    ];
    if x <= 2.0 {
        let y = x * x / 4.0;
        (x / 2.0).ln() * besseli1(x) + (1.0 / x) * poly(y, &P)
    } else {
        let y = 2.0 / x;
        (x.exp().recip() / x.sqrt()) * poly(y, &Q)
    }
}

/// Modified Bessel function Kn by upward recurrence (upstream `besselk`).
///
/// Upstream starts from `K0`, `K1` and applies
/// `K_{k+1} = K_{k-1} + k * (2/x) * K_k` for `k = 1..n-1`, so for `n <= 1`
/// the function returns `K1(x)` (reproduced exactly for fidelity).
pub fn besselk(n: i32, x: f64) -> f64 {
    let tox = 2.0 / x;
    let mut bkm = besselk0(x);
    let mut bk = besselk1(x);
    for j in 1..n {
        let bkp = bkm + (j as f64) * tox * bk;
        bkm = bk;
        bk = bkp;
    }
    bk
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Upstream fixture (reference/fixtures/fortran/test_bessel.txt) generated
    /// by `scripts/build_fortran_fixtures.sh` from jadexter/grtrans @ c76cb11.
    fn fixture() -> Vec<(f64, f64, f64, f64, f64)> {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../reference/fixtures/fortran/test_bessel.txt"
        );
        let text = std::fs::read_to_string(path)
            .expect("bessel fixture missing; run scripts/build_fortran_fixtures.sh");
        let mut out = Vec::new();
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let vals: Vec<f64> = line
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            if vals.len() != 5 {
                continue; // kn rows handled separately
            }
            out.push((vals[0], vals[1], vals[2], vals[3], vals[4]));
        }
        out
    }

    #[test]
    fn matches_upstream_fixture() {
        for (x, i0, i1, k0, k1) in fixture() {
            let (fi0, fi1, fk0, fk1) = (besseli0(x), besseli1(x), besselk0(x), besselk1(x));
            let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
            assert!(rel(fi0, i0) < 1e-14, "I0({x}): {fi0} vs {i0}");
            assert!(rel(fi1, i1) < 1e-14, "I1({x}): {fi1} vs {i1}");
            assert!(rel(fk0, k0) < 1e-14, "K0({x}): {fk0} vs {k0}");
            assert!(rel(fk1, k1) < 1e-14, "K1({x}): {fk1} vs {k1}");
        }
    }

    #[test]
    fn kn_matches_upstream_fixture() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../reference/fixtures/fortran/test_bessel.txt"
        );
        let text = std::fs::read_to_string(path).unwrap();
        let mut in_kn = false;
        for line in text.lines() {
            if line.starts_with("# kn") {
                in_kn = true;
                continue;
            }
            if !in_kn || line.trim().is_empty() {
                continue;
            }
            let vals: Vec<f64> = line
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            let (x, k2, k3, k5) = (vals[0], vals[1], vals[2], vals[3]);
            let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
            assert!(rel(besselk(2, x), k2) < 1e-14, "K2({x})");
            assert!(rel(besselk(3, x), k3) < 1e-14, "K3({x})");
            assert!(rel(besselk(5, x), k5) < 1e-14, "K5({x})");
        }
    }
}
