//! Polarized radiative transfer integration.
//!
//! Direct translation of `radtrans_integrate.f90` (upstream GRTRANS):
//! the Del Zanna & Bucciantini (2002) scheme (`delo`), the formal solution
//! (`formal`), the intensity-only quadrature, and the shared 4x4 matrix
//! helpers. The LSODA path is provided by [`lsoda`].
// Lints deliberately allowed: the port preserves upstream Fortran
// semantics and constants exactly.
// * `approx_constant`: upstream uses rounded literals (e.g. the
//   Euler-Mascheroni value -0.57721566 in bessel.f90) that must not be
//   replaced by the standard constants.
// * `neg_cmp_op_on_partial_ord`: `!(x > 0.0)` reproduces Fortran
//   `merge(a, b, x > 0)` semantics (NaN takes the false branch), which
//   `x <= 0.0` would not.
// * `too_many_arguments`: upstream routine signatures are preserved.
#![allow(clippy::approx_constant)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)] // index loops mirror Fortran 1-based indexing
#![allow(clippy::chunks_exact_to_as_chunks)]
#![allow(clippy::type_complexity)]
#![allow(clippy::needless_update)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_swap)]
#![allow(clippy::let_and_return)]
#![allow(clippy::items_after_test_module)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::assign_op_pattern)]
#![allow(unused_assignments)] // Fortran re-assignment patterns are preserved

pub mod lsoda;

/// Integration method (upstream `iflag`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// upstream iflag=0 (linear Stokes LSODA)
    Lsoda,
    /// upstream iflag=1
    Delo,
    /// upstream iflag=2
    Formal,
    /// upstream iflag=3 (experimental spherical Stokes LSODA)
    LsodaSph,
}

impl Method {
    pub fn from_name(iname: &str) -> Self {
        match iname {
            "lsoda" => Method::Lsoda,
            "delo" => Method::Delo,
            "formal" => Method::Formal,
            "lsodasph" => Method::LsodaSph,
            _ => panic!("iname not recognized: {iname}"),
        }
    }
}

/// 4x4 identity.
pub fn identity4() -> [[f64; 4]; 4] {
    let mut m = [[0.0f64; 4]; 4];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    m
}

/// Upstream `imatrix_4_single`: 4x4 inverse via the cofactor expansion.
pub fn imatrix_4(a: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let a11 = a[0][0];
    let a12 = a[0][1];
    let a13 = a[0][2];
    let a14 = a[0][3];
    let a21 = a[1][0];
    let a22 = a[1][1];
    let a23 = a[1][2];
    let a24 = a[1][3];
    let a31 = a[2][0];
    let a32 = a[2][1];
    let a33 = a[2][2];
    let a34 = a[2][3];
    let a41 = a[3][0];
    let a42 = a[3][1];
    let a43 = a[3][2];
    let a44 = a[3][3];
    let mut b = [[0.0f64; 4]; 4];
    b[0][0] = a22 * a33 * a44 + a23 * a34 * a42 + a24 * a32 * a43
        - a22 * a34 * a43
        - a23 * a32 * a44
        - a24 * a33 * a42;
    b[0][1] = a12 * a34 * a43 + a13 * a32 * a44 + a14 * a33 * a42
        - a12 * a33 * a44
        - a13 * a34 * a42
        - a14 * a32 * a43;
    b[0][2] = a12 * a23 * a44 + a13 * a24 * a42 + a14 * a22 * a43
        - a12 * a24 * a43
        - a13 * a22 * a44
        - a14 * a23 * a42;
    b[0][3] = a12 * a24 * a33 + a13 * a22 * a34 + a14 * a23 * a32
        - a12 * a23 * a34
        - a13 * a24 * a32
        - a14 * a22 * a33;
    b[1][0] = a21 * a34 * a43 + a23 * a31 * a44 + a24 * a33 * a41
        - a21 * a33 * a44
        - a23 * a34 * a41
        - a24 * a31 * a43;
    b[1][1] = a11 * a33 * a44 + a13 * a34 * a41 + a14 * a31 * a43
        - a11 * a34 * a43
        - a13 * a31 * a44
        - a14 * a33 * a41;
    b[1][2] = a11 * a24 * a43 + a13 * a21 * a44 + a14 * a23 * a41
        - a11 * a23 * a44
        - a13 * a24 * a41
        - a14 * a21 * a43;
    b[1][3] = a11 * a23 * a34 + a13 * a24 * a31 + a14 * a21 * a33
        - a11 * a24 * a33
        - a13 * a21 * a34
        - a14 * a23 * a31;
    b[2][0] = a21 * a32 * a44 + a22 * a34 * a41 + a24 * a31 * a42
        - a21 * a34 * a42
        - a22 * a31 * a44
        - a24 * a32 * a41;
    b[2][1] = a11 * a34 * a42 + a12 * a31 * a44 + a14 * a32 * a41
        - a11 * a32 * a44
        - a12 * a34 * a41
        - a14 * a31 * a42;
    b[2][2] = a11 * a22 * a44 + a12 * a24 * a41 + a14 * a21 * a42
        - a11 * a24 * a42
        - a12 * a21 * a44
        - a14 * a22 * a41;
    b[2][3] = a11 * a24 * a32 + a12 * a21 * a34 + a14 * a22 * a31
        - a11 * a22 * a34
        - a12 * a24 * a31
        - a14 * a21 * a32;
    b[3][0] = a21 * a33 * a42 + a22 * a31 * a43 + a23 * a32 * a41
        - a21 * a32 * a43
        - a22 * a33 * a41
        - a23 * a31 * a42;
    b[3][1] = a11 * a32 * a43 + a12 * a33 * a41 + a13 * a31 * a42
        - a11 * a33 * a42
        - a12 * a31 * a43
        - a13 * a32 * a41;
    b[3][2] = a11 * a23 * a42 + a12 * a21 * a43 + a13 * a22 * a41
        - a11 * a22 * a43
        - a12 * a23 * a41
        - a13 * a21 * a42;
    b[3][3] = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32
        - a11 * a23 * a32
        - a12 * a21 * a33
        - a13 * a22 * a31;
    let det_a = a[0][0] * b[0][0] + a[1][0] * b[0][1] + a[2][0] * b[0][2] + a[3][0] * b[0][3];
    for row in b.iter_mut() {
        for v in row.iter_mut() {
            *v /= det_a;
        }
    }
    b
}

/// Upstream `opacity_matrix`.
pub fn opacity_matrix(a: &[f64], p: &[f64]) -> [[f64; 4]; 4] {
    let mut k = [[0.0f64; 4]; 4];
    for i in 0..4 {
        k[i][i] = a[0];
    }
    k[1][0] = a[1];
    k[2][0] = a[2];
    k[3][0] = a[3];
    k[0][1] = a[1];
    k[0][2] = a[2];
    k[0][3] = a[3];
    k[2][1] = -p[2];
    k[3][1] = p[1];
    k[1][2] = p[2];
    k[1][3] = -p[1];
    k[3][2] = -p[0];
    k[2][3] = p[0];
    k
}

fn matvec(m: &[[f64; 4]; 4], v: &[f64; 4]) -> [f64; 4] {
    let mut out = [0.0f64; 4];
    for i in 0..4 {
        out[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2] + m[i][3] * v[3];
    }
    out
}

fn matmul(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut s = 0.0;
            for k in 0..4 {
                s += a[i][k] * b[k][j];
            }
            out[i][j] = s;
        }
    }
    out
}

/// Del Zanna & Bucciantini (2002) integration scheme
/// (upstream `radtrans_integrate_delo`).
///
/// * `x`: affine parameter samples (npts)
/// * `tau`: optical depth samples (npts)
/// * `j`: emission (npts x 4)
/// * `a`: absorption terms (npts x 4: alpha_I,Q,U,V)
/// * `rho`: Faraday terms (npts x 3: rho_Q,U,V)
/// * `thin`: threshold between the thick and thin branches (upstream `thin`)
///
/// Returns the intensity (4 x npts).
#[allow(clippy::too_many_arguments)]
pub fn radtrans_integrate_delo(
    x: &[f64],
    tau: &[f64],
    j: &[[f64; 4]],
    a: &[[f64; 4]],
    rho: &[[f64; 3]],
    thin: f64,
) -> Vec<[f64; 4]> {
    let npts = x.len();
    let identity = identity4();
    let mut intensity = vec![[0.0f64; 4]; npts];
    let delta: Vec<f64> = (1..npts).map(|k| tau[k] - tau[k - 1]).collect();
    let karr: Vec<[[f64; 4]; 4]> = (0..npts).map(|k| opacity_matrix(&a[k], &rho[k])).collect();
    let dx: Vec<f64> = (1..npts).map(|k| x[k - 1] - x[k]).collect();
    let mut iprev = [0.0f64; 4];
    intensity[0] = iprev;
    for k in (0..npts - 1).rev() {
        let k0 = karr[k];
        let k1 = karr[k + 1];
        let (pt, qt);
        if delta[k] > thin {
            let e = (-delta[k]).exp();
            let f = 1.0 - e;
            let g = (1.0 - (1.0 + delta[k]) * e) / delta[k];
            let sp: [f64; 4] = [
                j[k][0] / a[k][0],
                j[k][1] / a[k][0],
                j[k][2] / a[k][0],
                j[k][3] / a[k][0],
            ];
            let sp1: [f64; 4] = [
                j[k + 1][0] / a[k + 1][0],
                j[k + 1][1] / a[k + 1][0],
                j[k + 1][2] / a[k + 1][0],
                j[k + 1][3] / a[k + 1][0],
            ];
            let mut kp = k0;
            let mut kp1 = k1;
            for i in 0..4 {
                for q in 0..4 {
                    kp[i][q] /= a[k][0];
                    kp1[i][q] /= a[k + 1][0];
                }
                kp[i][i] -= 1.0;
                kp1[i][i] -= 1.0;
            }
            let mut matrix = identity;
            for i in 0..4 {
                for q in 0..4 {
                    matrix[i][q] += (f - g) * kp[i][q];
                }
            }
            let imatrix = imatrix_4(&matrix);
            let mut pv = [0.0f64; 4];
            for i in 0..4 {
                pv[i] = (f - g) * sp[i] + g * sp1[i];
            }
            pt = matvec(&imatrix, &pv);
            let mut qm = identity;
            for i in 0..4 {
                for q in 0..4 {
                    qm[i][q] = e * identity[i][q] - g * kp1[i][q];
                }
            }
            qt = matmul(&imatrix, &qm);
        } else {
            let ki = a[k][0];
            let d = delta[k];
            let mut matrix = identity;
            for i in 0..4 {
                for q in 0..4 {
                    matrix[i][q] = (1.0 - d / 2.0 + d * d / 6.0) * identity[i][q]
                        + (0.5 * dx[k] - dx[k] * dx[k] / 6.0 * ki) * k0[i][q];
                }
            }
            let imatrix = imatrix_4(&matrix);
            let mut pv = [0.0f64; 4];
            for i in 0..4 {
                pv[i] = (0.5 * dx[k] * j[k][i] - dx[k] * dx[k] / 6.0 * ki * j[k][i])
                    + (0.5 * j[k + 1][i] * dx[k] - dx[k] * dx[k] / 3.0 * ki * j[k + 1][i]);
            }
            pt = matvec(&imatrix, &pv);
            let mut qm = identity;
            for i in 0..4 {
                for q in 0..4 {
                    qm[i][q] = identity[i][q]
                        * (1.0 - 0.5 * dx[k] * ki + dx[k] * dx[k] / 6.0 * ki * ki)
                        - (0.5 * dx[k] - dx[k] * dx[k] / 3.0) * k1[i][q];
                }
            }
            qt = matmul(&imatrix, &qm);
        }
        let mut newi = [0.0f64; 4];
        let tmp = matvec(&qt, &iprev);
        for i in 0..4 {
            newi[i] = pt[i] + tmp[i];
        }
        intensity[npts - k - 1] = newi;
        iprev = newi;
    }
    intensity
}

/// Upstream `calc_O` (formal-solution matrix exponential).
#[allow(clippy::too_many_arguments)]
pub fn calc_o(
    a: &[f64; 4],
    rho: &[f64; 3],
    dx: f64,
    identity: &[[f64; 4]; 4],
    m1: &mut [[f64; 4]; 4],
    m2: &mut [[f64; 4]; 4],
    m3: &mut [[f64; 4]; 4],
    m4: &mut [[f64; 4]; 4],
) -> [[f64; 4]; 4] {
    let onopol = (-a[0] * dx).exp();
    let aq = a[1];
    let au = a[2];
    let av = a[3];
    let rhoq = rho[0];
    let rhou = rho[1];
    let rhov = rho[2];
    let a2 = aq * aq + au * au + av * av;
    let p2 = rhoq * rhoq + rhou * rhou + rhov * rhov;
    if a2 == 0.0 && p2 == 0.0 {
        let mut o = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for q in 0..4 {
                o[i][q] = identity[i][q] * onopol;
            }
        }
        return o;
    }
    let ap = aq * rhoq + au * rhou + av * rhov;
    let disc = ((a2 - p2) * (a2 - p2) / 4.0 + ap * ap).sqrt();
    let lam1 = (disc + (a2 - p2) / 2.0).sqrt();
    let lam2 = (disc - (a2 - p2) / 2.0).sqrt();
    let theta = lam1 * lam1 + lam2 * lam2;
    let sig = if ap >= 0.0 { 1.0 } else { -1.0 };
    let m1v = *identity;
    let mut m2v = [[0.0f64; 4]; 4];
    m2v[1][0] = lam2 * aq - sig * lam1 * rhoq;
    m2v[2][0] = lam2 * au - sig * lam1 * rhou;
    m2v[3][0] = lam2 * av - sig * lam1 * rhov;
    m2v[0][1] = lam2 * aq - sig * lam1 * rhoq;
    m2v[2][1] = -sig * lam1 * av - lam2 * rhov;
    m2v[3][1] = sig * lam1 * au + lam2 * rhou;
    m2v[0][2] = lam2 * au - sig * lam1 * rhou;
    m2v[1][2] = sig * lam1 * av + lam2 * rhov;
    m2v[3][2] = -sig * lam1 * aq - lam2 * rhoq;
    m2v[0][3] = lam2 * av - sig * lam1 * rhov;
    m2v[1][3] = -sig * lam1 * au - lam2 * rhou;
    m2v[2][3] = sig * lam1 * aq + lam2 * rhoq;
    for row in m2v.iter_mut() {
        for v in row.iter_mut() {
            *v /= theta;
        }
    }
    let mut m3v = [[0.0f64; 4]; 4];
    m3v[1][0] = lam1 * aq + sig * lam2 * rhoq;
    m3v[2][0] = lam1 * au + sig * lam2 * rhou;
    m3v[3][0] = lam1 * av + sig * lam2 * rhov;
    m3v[0][1] = lam1 * aq + sig * lam2 * rhoq;
    m3v[2][1] = sig * lam2 * av - lam1 * rhov;
    m3v[3][1] = -sig * lam2 * au + lam1 * rhou;
    m3v[0][2] = lam1 * au + sig * lam2 * rhou;
    m3v[1][2] = -sig * lam2 * av + lam1 * rhov;
    m3v[3][2] = sig * lam2 * aq - lam1 * rhoq;
    m3v[0][3] = lam1 * av + sig * lam2 * rhov;
    m3v[1][3] = sig * lam2 * au - lam1 * rhou;
    m3v[2][3] = -sig * lam2 * aq + lam1 * rhoq;
    for row in m3v.iter_mut() {
        for v in row.iter_mut() {
            *v /= theta;
        }
    }
    let mut m4v = [[0.0f64; 4]; 4];
    m4v[0][0] = (a2 + p2) / 2.0;
    m4v[1][0] = au * rhov - av * rhou;
    m4v[2][0] = av * rhoq - aq * rhov;
    m4v[3][0] = aq * rhou - au * rhoq;
    m4v[0][1] = av * rhou - au * rhov;
    m4v[1][1] = aq * aq + rhoq * rhoq - (a2 + p2) / 2.0;
    m4v[2][1] = aq * au + rhoq * rhou;
    m4v[3][1] = av * aq + rhov * rhoq;
    m4v[0][2] = aq * rhov - av * rhoq;
    m4v[1][2] = aq * au + rhoq * rhou;
    m4v[2][2] = au * au + rhou * rhou - (a2 + p2) / 2.0;
    m4v[3][2] = au * av + rhou * rhov;
    m4v[0][3] = au * rhoq - aq * rhou;
    m4v[1][3] = av * aq + rhov * rhoq;
    m4v[2][3] = au * av + rhou * rhov;
    m4v[3][3] = av * av + rhov * rhov - (a2 + p2) / 2.0;
    for row in m4v.iter_mut() {
        for v in row.iter_mut() {
            *v *= 2.0 / theta;
        }
    }
    *m1 = m1v;
    *m2 = m2v;
    *m3 = m3v;
    *m4 = m4v;
    let mut o = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for q in 0..4 {
            o[i][q] = onopol
                * (0.5 * ((lam1 * dx).cosh() + (lam2 * dx).cos()) * m1v[i][q]
                    - (lam2 * dx).sin() * m2v[i][q]
                    - (lam1 * dx).sinh() * m3v[i][q]
                    + 0.5 * ((lam1 * dx).cosh() - (lam2 * dx).cos()) * m4v[i][q]);
        }
    }
    o
}

/// Upstream `radtrans_integrate_formal`.
pub fn radtrans_integrate_formal(
    x: &[f64],
    j: &[[f64; 4]],
    a: &[[f64; 4]],
    rho: &[[f64; 3]],
) -> Vec<[f64; 4]> {
    let npts = x.len();
    let identity = identity4();
    let mut intensity = vec![[0.0f64; 4]; npts];
    let dx: Vec<f64> = (1..npts).map(|k| x[k - 1] - x[k]).collect();
    let mut iprev = [0.0f64; 4];
    intensity[0] = iprev;
    let mut m1 = identity;
    let mut m2 = identity;
    let mut m3 = identity;
    let mut m4 = identity;
    for k in (0..npts - 1).rev() {
        let o = calc_o(
            &a[k], &rho[k], dx[k], &identity, &mut m1, &mut m2, &mut m3, &mut m4,
        );
        let mut v = [0.0f64; 4];
        for i in 0..4 {
            v[i] = j[k][i] * dx[k] + iprev[i];
        }
        let inew = matvec(&o, &v);
        iprev = inew;
        intensity[npts - k - 1] = inew;
    }
    intensity
}

/// Upstream `radtrans_integrate_quadrature` (intensity only, I0 = 0).
pub fn radtrans_integrate_quadrature(
    s: &[f64],
    j: &[f64],
    _kcoef: &[f64],
    tau: &[f64],
) -> Vec<f64> {
    let n = s.len();
    let integrand: Vec<f64> = (0..n).map(|i| j[i] * (-tau[i]).exp()).collect();
    // intensity(1,:) = -tsum(s, j*exp(-tau))
    let ts = grtrans_core::math::tsum(s, &integrand);
    ts.iter().map(|v| -v).collect()
}

/// Upstream `calc_opt_depth`: cumulative trapezoidal optical depth.
pub fn calc_opt_depth(s: &[f64], alpha_i: &[f64]) -> Vec<f64> {
    let abs: Vec<f64> = alpha_i.iter().map(|v| v.abs()).collect();
    grtrans_core::math::tsum(s, &abs)
}

/// Driver-facing `integrate` dispatch for the intensity-only case
/// (`nequations==1`): upstream `lsoda` solves the scalar ODE, while
/// `delo`/`formal` use the trapezoidal quadrature.
pub fn integrate_scalar(
    method: Method,
    s: &[f64],
    j: &[f64],
    k: &[f64],
    tau: &[f64],
    hmax: f64,
    oatol: f64,
    ortol: f64,
) -> (Vec<f64>, usize) {
    match method {
        Method::Lsoda => lsoda::integrate_lsoda_scalar(s, j, k, tau, oatol, ortol, hmax),
        _ => {
            let n = s.len();
            (radtrans_integrate_quadrature(s, j, k, tau), n)
        }
    }
}

/// Driver-facing `integrate` dispatch for the 4-Stokes case.
#[allow(clippy::too_many_arguments)]
pub fn integrate(
    method: Method,
    s: &[f64],
    j: &[[f64; 4]],
    kcoef: &[[f64; 7]],
    tau: &[f64],
    thin: f64,
    hmax: f64,
    oatol: f64,
    ortol: f64,
) -> (Vec<[f64; 4]>, usize) {
    let a: Vec<[f64; 4]> = kcoef.iter().map(|k| [k[0], k[1], k[2], k[3]]).collect();
    let rho: Vec<[f64; 3]> = kcoef.iter().map(|k| [k[4], k[5], k[6]]).collect();
    match method {
        Method::Delo => {
            let n = s.len();
            (radtrans_integrate_delo(s, tau, j, &a, &rho, thin), n)
        }
        Method::Formal => {
            let n = s.len();
            (radtrans_integrate_formal(s, j, &a, &rho), n)
        }
        Method::Lsoda => lsoda::integrate_lsoda(s, j, &a, &rho, tau, oatol, ortol, hmax),
        Method::LsodaSph => panic!("lsodasph is experimental upstream and not ported"),
    }
}
