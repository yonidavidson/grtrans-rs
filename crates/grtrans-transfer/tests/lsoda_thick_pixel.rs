//! Optically thick SARIAF pixel: our stiff-stable integrator must match an
//! independent high-accuracy RK4 solution of the upstream shifted ODE, NOT
//! upstream LSODA (which is itself ~20% off on this pixel).
//!
//! Fixture: upstream-dumped j/K/s/tau and LSODA intensity for one pixel with
//! K up to ~5e2 and tau up to ~1e3. See docs/VALIDATION_PLAN.md.
#![allow(clippy::chunks_exact_to_as_chunks)]

use grtrans_transfer::lsoda::integrate_lsoda;
use grtrans_transfer::opacity_matrix;

fn read(name: &str) -> Vec<f64> {
    let p = format!(
        "{}/../../reference/fixtures/sariaf/pixel/{name}.f64.bin",
        env!("CARGO_MANIFEST_DIR")
    );
    let b = std::fs::read(&p).unwrap_or_else(|_| panic!("missing {p}"));
    b.chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn shifted<const N: usize>(s: &[f64], vals: &[[f64; N]], lam: f64) -> [f64; N] {
    let n = s.len();
    let mut lo = 0usize;
    let mut hi = n - 1;
    if lam <= s[n - 1] {
        lo = 0;
    } else if lam >= s[0] {
        lo = n - 2;
    } else {
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if s[n - 1 - mid] <= lam {
                lo = mid;
            } else {
                hi = mid;
            }
        }
    }
    let l = lo;
    let ss_l = s[n - 1 - l];
    let ss_l1 = s[n - 1 - (l + 1)];
    let w = (lam - ss_l) / (ss_l1 - ss_l);
    let indx0 = n - 1 - l;
    let uindx0 = (indx0 + 1).min(n - 1);
    let mut out = [0.0f64; N];
    for c in 0..N {
        out[c] = (1.0 - w) * vals[indx0][c] + w * vals[uindx0][c];
    }
    out
}

fn rhs(
    s: &[f64],
    j: &[[f64; 4]],
    a: &[[f64; 4]],
    rho: &[[f64; 3]],
    lam: f64,
    y: &[f64; 4],
) -> [f64; 4] {
    let jv = shifted(s, j, lam);
    let av = shifted(s, a, lam);
    let rv = shifted(s, rho, lam);
    let k = opacity_matrix(&av, &rv);
    let mut out = [0.0f64; 4];
    for i in 0..4 {
        out[i] = jv[i] - (k[i][0] * y[0] + k[i][1] * y[1] + k[i][2] * y[2] + k[i][3] * y[3]);
    }
    out
}

fn rk4_exact(s: &[f64], j: &[[f64; 4]], a: &[[f64; 4]], rho: &[[f64; 3]], iobs: usize) -> f64 {
    let lam_start = s[iobs];
    let lam_end = s[0];
    let nsteps = 400000usize;
    let h = (lam_end - lam_start) / nsteps as f64;
    let mut y = [0.0f64; 4];
    let mut lam = lam_start;
    let add = |u: &[f64; 4], v: &[f64; 4], f: f64| -> [f64; 4] {
        [
            u[0] + f * v[0],
            u[1] + f * v[1],
            u[2] + f * v[2],
            u[3] + f * v[3],
        ]
    };
    for _ in 0..nsteps {
        let k1 = rhs(s, j, a, rho, lam, &y);
        let y2 = add(&y, &k1, h / 2.0);
        let k2 = rhs(s, j, a, rho, lam + h / 2.0, &y2);
        let y3 = add(&y, &k2, h / 2.0);
        let k3 = rhs(s, j, a, rho, lam + h / 2.0, &y3);
        let y4 = add(&y, &k3, h);
        let k4 = rhs(s, j, a, rho, lam + h, &y4);
        for i in 0..4 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        lam += h;
    }
    y[0]
}

#[test]
fn thick_pixel_matches_exact_not_lsoda_error() {
    let s = read("s");
    let jb = read("j");
    let kb = read("K");
    let tau = read("tau");
    let i1up = read("I1");
    let n = s.len();
    let mut j = vec![[0.0f64; 4]; n];
    let mut a = vec![[0.0f64; 4]; n];
    let mut rho = vec![[0.0f64; 3]; n];
    for i in 0..n {
        for c in 0..4 {
            j[i][c] = jb[c * n + i];
            a[i][c] = kb[c * n + i];
        }
        for c in 0..3 {
            rho[i][c] = kb[(4 + c) * n + i];
        }
    }
    let (ours, nptsout) = integrate_lsoda(&s, &j, &a, &rho, &tau, 1e-8, 1e-6, 10.0);
    let iobs = nptsout - 1;
    let exact = rk4_exact(&s, &j, &a, &rho, iobs);
    let rel = (ours[iobs][0] - exact).abs() / exact.abs();
    println!(
        "thick pixel: ours={:.10e} exact={:.10e} upstream={:.10e} rel(ours,exact)={:.3e} upstream_dev={:.4}",
        ours[iobs][0],
        exact,
        i1up[iobs],
        rel,
        ours[iobs][0] / i1up[iobs] - 1.0
    );
    assert!(
        rel < 1e-4,
        "integrator deviates from exact solution: {rel:e}"
    );
}
