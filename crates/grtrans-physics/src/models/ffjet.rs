//! Semi-analytic jet model of Broderick & Loeb (2009).
//!
//! Direct translation of `fluid_model_ffjet.f90` (upstream GRTRANS). The
//! model reads a Fortran-unformatted binary grid (nx = nx1*nx2, r fastest)
//! and bilinearly interpolates rho, internal energy, field components and
//! velocities; the four-velocity is reconstructed with the (single
//! precision, as resolved by gfortran for the mixed-kind call) inverse LNRF
//! transform.
//!
//! Note: gfortran resolves the upstream mixed-kind `kerr_metric(zr,
//! real(th), a)` call to the single-precision overload, so the metric used
//! for bmag and the four-velocity is computed in f32 and widened; the Rust
//! port reproduces that.

use grtrans_core::four_vector::FourVector;
use grtrans_core::interpolate::bilin_f32;
use grtrans_core::kerr::{blmetric_cov, blmetric_cov_f32, lnrf_frame_inv_f32};
use std::path::Path;

/// Upstream module state (grids and fluid arrays).
#[derive(Clone, Debug, Default)]
pub struct FfjetData {
    pub nx1: usize,
    pub nx2: usize,
    pub rc_arr: Vec<f32>,
    pub thc_arr: Vec<f32>,
    pub rho_arr: Vec<f32>,
    pub p_arr: Vec<f32>,
    pub u0_arr: Vec<f32>,
    pub vr_arr: Vec<f32>,
    pub vth_arr: Vec<f32>,
    pub vph_arr: Vec<f32>,
    pub b0_arr: Vec<f32>,
    pub br_arr: Vec<f32>,
    pub bth_arr: Vec<f32>,
    pub bph_arr: Vec<f32>,
}

/// Read one Fortran unformatted record (4-byte length markers) as f32.
fn read_record_f32(data: &[u8], off: &mut usize) -> Vec<f32> {
    let m = i32::from_le_bytes(data[*off..*off + 4].try_into().unwrap()) as usize;
    *off += 4;
    let payload = &data[*off..*off + m];
    let vals: Vec<f32> = payload
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
        .collect();
    *off += m + 4;
    vals
}

/// Upstream `initialize_ffjet_model`: read the binary dump.
pub fn initialize_ffjet_model(path: &Path) -> FfjetData {
    let data = std::fs::read(path)
        .unwrap_or_else(|e| panic!("cannot read FFJET dump {}: {e}", path.display()));
    let mut off = 0usize;
    // record 0: aa (f32), nx (i32)
    let m0 = i32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as usize;
    off += 4;
    let aa = f32::from_le_bytes(data[off..off + 4].try_into().unwrap());
    let nx = i32::from_le_bytes(data[off + 4..off + 8].try_into().unwrap()) as usize;
    off += m0 + 4;
    let _ = aa;
    let nx1 = (nx as f64).sqrt() as usize;
    assert_eq!(nx1 * nx1, nx, "FFJET grid is not square");
    let rec1 = read_record_f32(&data, &mut off);
    let rec2 = read_record_f32(&data, &mut off);
    let rec3 = read_record_f32(&data, &mut off);
    let rec4 = read_record_f32(&data, &mut off);
    assert_eq!(rec1.len(), 3 * nx);
    assert_eq!(rec2.len(), 5 * nx);
    assert_eq!(rec3.len(), 4 * nx);
    assert_eq!(rec4.len(), 2 * nx);
    FfjetData {
        nx1,
        nx2: nx1,
        rc_arr: rec1[0..nx].to_vec(),
        thc_arr: rec1[nx..2 * nx].to_vec(),
        rho_arr: rec1[2 * nx..3 * nx].to_vec(),
        // rec2: b (dummy), b0, br, bth, bph
        b0_arr: rec2[nx..2 * nx].to_vec(),
        br_arr: rec2[2 * nx..3 * nx].to_vec(),
        bth_arr: rec2[3 * nx..4 * nx].to_vec(),
        bph_arr: rec2[4 * nx..5 * nx].to_vec(),
        u0_arr: rec3[0..nx].to_vec(),
        vr_arr: rec3[nx..2 * nx].to_vec(),
        vth_arr: rec3[2 * nx..3 * nx].to_vec(),
        vph_arr: rec3[3 * nx..4 * nx].to_vec(),
        p_arr: vec![0.0; nx],
    }
}

/// Upstream `ffjet_vals`: interpolate the jet solution at each point.
/// Returns `(rho, p, bmag, u, b)` per point.
pub fn ffjet_vals(
    data: &FfjetData,
    x0: &[FourVector],
    a: f32,
) -> (
    Vec<f32>,
    Vec<f32>,
    Vec<f32>,
    Vec<FourVector>,
    Vec<FourVector>,
) {
    let nx1 = data.nx1;
    let nx2 = data.nx2;
    let npts = x0.len();
    let umax = (nx2 - 1) as i32;
    let uniqr: Vec<f32> = data.rc_arr[0..nx1].to_vec();
    let uniqx2: Vec<f32> = (0..nx1).map(|i| data.thc_arr[i * nx1]).collect();
    let uniqth = uniqx2.clone();
    let uniqx1: Vec<f32> = uniqr.iter().map(|v| v.ln()).collect();

    let mut rho = vec![0.0f32; npts];
    let mut p = vec![0.0f32; npts];
    let mut bmag = vec![0.0f32; npts];
    let mut u = vec![FourVector::flat([0.0; 4]); npts];
    let mut b = vec![FourVector::flat([0.0; 4]); npts];

    for i in 0..npts {
        // upstream: zm = cos(theta) stored in a default real; x2 = acos(|zm|)
        // in single precision; zr = x0%data(2) stored in a default real
        let zm = (x0[i].data[2].cos()) as f32;
        let x2 = zm.abs().acos();
        let zr = x0[i].data[1] as f32;
        let x1 = zr.ln();
        let mut lx1 =
            ((x1 - uniqx1[0]) / (uniqx1[nx1 - 1] - uniqx1[0]) * (nx1 as f32 - 1.0)) as i32 + 1;
        let mut lx2 =
            ((x2 - uniqx2[0]) / (uniqx2[nx2 - 1] - uniqx2[0]) * (nx2 as f32 - 1.0)) as i32 + 1;
        // upstream clamps lx2 only; lx1 is clamped here for memory safety
        // (out-of-range rays read arbitrary in-bounds memory upstream)
        lx1 = lx1.clamp(1, nx1 as i32 - 1);
        lx2 = lx2.clamp(1, umax);
        let ux1 = lx1 + 1;
        let ux2 = lx2 + 1;
        let td = (x2 - uniqth[(lx2 - 1) as usize])
            / (uniqth[(ux2 - 1) as usize] - uniqth[(lx2 - 1) as usize]);
        let rd = (zr - uniqr[(lx1 - 1) as usize])
            / (uniqr[(ux1 - 1) as usize] - uniqr[(lx1 - 1) as usize]);
        let x1l = lx1 - 1;
        let x1u = ux1 - 1;
        let x2l = (lx2 - 1) * nx1 as i32;
        let x2u = (ux2 - 1) * nx1 as i32;
        let idx = [
            (x1l + x2l) as usize,
            (x1l + x2u) as usize,
            (x1u + x2l) as usize,
            (x1u + x2u) as usize,
        ];
        let get4 =
            |arr: &[f32]| -> [f32; 4] { [arr[idx[0]], arr[idx[1]], arr[idx[2]], arr[idx[3]]] };
        let bil = |v: [f32; 4]| -> f32 { bilin_f32(v[0], v[1], v[2], v[3], rd, td) };
        let inside = x1 > uniqx1[0];
        let rho_v = bil(get4(&data.rho_arr));
        rho[i] = if inside { rho_v } else { 0.0 };
        let p_v = bil(get4(&data.p_arr));
        p[i] = if inside { p_v } else { 1.0 };
        let vrl0 = if inside { bil(get4(&data.vr_arr)) } else { 0.0 };
        let vtl0 = if inside {
            bil(get4(&data.vth_arr))
        } else {
            1.0
        };
        let vpl0 = if inside {
            bil(get4(&data.vph_arr))
        } else {
            0.0
        };
        for (c, arr) in [
            (0usize, &data.b0_arr),
            (1, &data.br_arr),
            (2, &data.bth_arr),
            (3, &data.bph_arr),
        ] {
            let v = bil(get4(arr));
            b[i].data[c] = if inside { v as f64 } else { 1.0 };
        }
        let u0_v = bil(get4(&data.u0_arr));
        u[i].data[0] = if inside { u0_v as f64 } else { 1.0 };
        // metric for bmag: upstream call kerr_metric(zr, real(th), a) with
        // all-default-real arguments resolves to the single-precision
        // overload; values are widened to f64 for the four-vector storage
        let metric = blmetric_cov_f32(zr, x0[i].data[2] as f32, a);
        b[i].assign_metric(metric);
        let bdotb = b[i].dot(&b[i]);
        bmag[i] = bdotb.max(0.0).sqrt() as f32;
        // four-velocity: inverse LNRF transform (single precision)
        let (vr0, vth0, vph0) = lnrf_frame_inv_f32(vrl0, vtl0, vpl0, zr, a, x0[i].data[2] as f32);
        u[i].data[1] = u[i].data[0] * vr0 as f64;
        u[i].data[2] = u[i].data[0] * vth0 as f64;
        u[i].data[3] = u[i].data[0] * vph0 as f64;
        let metric64 = blmetric_cov(x0[i].data[1], x0[i].data[2], a as f64);
        u[i].assign_metric(metric64);
    }
    (rho, p, bmag, u, b)
}
