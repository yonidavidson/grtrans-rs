//! HARM GRMHD dump reader and interpolation.
//!
//! Direct translation of `fluid_model_harm.f90` (upstream GRTRANS) and the
//! `get_harm_fluidvars` / `convert_fluidvars_harm` / `scale_sim_units`
//! routines from `fluid.f90`. Only single-snapshot (`nt=1`) data is
//! supported, which is what the upstream regression case uses; multi-file
//! time interpolation is documented as not ported in
//! docs/PORTING_MATRIX.md.
//!
//! The data file is a formatted HARM dump (34 columns per line, `nx2` lines
//! per radial block). Columns (1-based): 1-4 = x1, x2, r, theta; 5 = rho;
//! 6 = pressure; 14-17 = u^t, u^r, u^theta, u^phi (MKS); 22-25 = b^t, b^r,
//! b^theta, b^phi (MKS); 34 = sqrt(-g).

use grtrans_core::constants::{C2, K, MP, MSUN, PI};
use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{
    blmetric_cov, blmetric_cov_f32, lnrf_frame_f32, lnrf_frame_inv_f32, uks2ubl,
};
use grtrans_core::math::zbrent;
use grtrans_core::constants::G as GC;
use std::path::Path;

/// Upstream module state for HARM.
#[derive(Clone, Debug, Default)]
pub struct HarmData {
    pub nx1: usize,
    pub nx2: usize,
    pub x1_arr: Vec<f32>,
    pub x2_arr: Vec<f32>,
    pub r_arr: Vec<f32>,
    pub th_arr: Vec<f32>,
    pub rho_arr: Vec<f32>,
    pub p_arr: Vec<f32>,
    /// BL u^t from the dump
    pub u0_arr: Vec<f32>,
    /// LNRF velocity components (computed on load, as upstream)
    pub vrl_arr: Vec<f32>,
    pub vtl_arr: Vec<f32>,
    pub vpl_arr: Vec<f32>,
    pub b0_arr: Vec<f32>,
    pub br_arr: Vec<f32>,
    pub bth_arr: Vec<f32>,
    pub bph_arr: Vec<f32>,
    pub h: f32,
    pub asim: f32,
    /// adiabatic index from the header (header(11)); the dump's pressure
    /// column is converted via p_arr = p*(gam-1) as upstream
    pub gam: f32,
}

/// Upstream `umksh2uks`: convert a MKS four-vector component set to KS
/// spherical components.
fn umksh2uks(
    fut: &FourVector,
    r: f32,
    x2: f32,
    hval: f32,
) -> FourVector {
    let dthdx2 = PI as f32 * (1.0 + (1.0 - hval) * (2.0 * PI as f32 * x2).cos());
    let mut out = *fut;
    out.data[1] = r as f64 * fut.data[1];
    out.data[2] = fut.data[2] * dthdx2 as f64;
    out
}

/// Upstream `transformbl2mksh`: BL -> MKS coordinate transform (x2 by root
/// find).
fn transformbl2mksh(r: f32, theta: f32, h: f32) -> (f32, f32) {
    let x1 = r.ln();
    let th = theta as f64;
    let hv = h as f64;
    let x2 = zbrent(
        |x| th - PI * x - 0.5 * (1.0 - hv) * (2.0 * PI * x).sin(),
        0.0,
        1.0,
        1e-6,
    );
    (x1, x2 as f32)
}

/// Read the HARM dump (single snapshot).
pub fn initialize_harm_model(path: &Path) -> HarmData {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read HARM dump {}: {e}", path.display()));
    let mut lines = text.lines();
    let header: Vec<f64> = lines
        .next()
        .expect("empty HARM dump")
        .split_whitespace()
        .take(26)
        .map(|t| t.replace(['d', 'D'], "e").parse().unwrap())
        .collect();
    assert_eq!(header.len(), 26, "HARM header must have 26 values");
    let nx1 = header[1] as usize;
    let nx2 = header[2] as usize;
    let asim = header[9] as f32;
    let h = header[24] as f32; // header(nhead-1), nhead=26
    let gam = header[10] as f32;
    let n = nx1 * nx2;
    let mut data = HarmData {
        nx1,
        nx2,
        h,
        asim,
        gam,
        ..Default::default()
    };
    let mut raw_u: Vec<[f32; 4]> = vec![[0.0; 4]; n];
    let mut raw_b: Vec<[f32; 4]> = vec![[0.0; 4]; n];
    let mut idx = 0usize;
    for line in lines {
        if idx >= n {
            break;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 34 {
            continue;
        }
        let v: Vec<f32> = toks[..34]
            .iter()
            .map(|t| t.replace(['d', 'D'], "e").parse().unwrap())
            .collect();
        data.x1_arr.push(v[0]);
        data.x2_arr.push(v[1]);
        data.r_arr.push(v[2]);
        data.th_arr.push(v[3]);
        data.rho_arr.push(v[4]);
        // the dump column is internal energy; upstream stores p*(gam-1)
        let p_raw: f64 = toks[5].replace(['d', 'D'], "e").parse().unwrap();
        data.p_arr.push((p_raw * (gam as f64 - 1.0)) as f32);
        raw_u[idx] = [v[13], v[14], v[15], v[16]];
        raw_b[idx] = [v[21], v[22], v[23], v[24]];
        idx += 1;
    }
    assert_eq!(idx, n, "HARM dump has {idx} of {n} rows");

    // transforms: MKS -> KS -> BL, then LNRF velocities
    data.u0_arr = vec![0.0; n];
    data.vrl_arr = vec![0.0; n];
    data.vtl_arr = vec![0.0; n];
    data.vpl_arr = vec![0.0; n];
    data.b0_arr = vec![0.0; n];
    data.br_arr = vec![0.0; n];
    data.bth_arr = vec![0.0; n];
    data.bph_arr = vec![0.0; n];
    for i in 0..n {
        let r = data.r_arr[i];
        let x2 = data.x2_arr[i];
        let th = data.th_arr[i];
        let fu = FourVector::flat(raw_u[i].map(f64::from));
        let uks = umksh2uks(&fu, r, x2, h);
        let ubl = uks2ubl(&uks, r as f64, asim as f64);
        let fb = FourVector::flat(raw_b[i].map(f64::from));
        let bks = umksh2uks(&fb, r, x2, h);
        let bbl = uks2ubl(&bks, r as f64, asim as f64);
        data.u0_arr[i] = ubl.data[0] as f32;
        let (vrl, vtl, vpl) = lnrf_frame_f32(
            (ubl.data[1] / ubl.data[0]) as f32,
            (ubl.data[2] / ubl.data[0]) as f32,
            (ubl.data[3] / ubl.data[0]) as f32,
            r,
            asim,
            th,
        );
        data.vrl_arr[i] = vrl;
        data.vtl_arr[i] = vtl;
        data.vpl_arr[i] = vpl;
        data.b0_arr[i] = bbl.data[0] as f32;
        data.br_arr[i] = bbl.data[1] as f32;
        data.bth_arr[i] = bbl.data[2] as f32;
        data.bph_arr[i] = bbl.data[3] as f32;
    }
    data
}

/// Upstream `harm_vals` for a single snapshot: interpolate the fluid
/// variables at each point. Returns `(rho, p, bmag, u, b)`.
pub fn harm_vals(
    data: &HarmData,
    x0: &[FourVector],
    a: f32,
) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<FourVector>, Vec<FourVector>) {
    let nx1 = data.nx1;
    let nx2 = data.nx2;
    let npts = x0.len();
    let uniqx2: Vec<f32> = data.x2_arr[0..nx2].to_vec();
    let uniqx1: Vec<f32> = (0..nx1).map(|i| data.x1_arr[i * nx2]).collect();
    let uniqr: Vec<f32> = uniqx1.iter().map(|v| v.exp()).collect();
    let uniqth: Vec<f32> = uniqx2
        .iter()
        .map(|v| PI as f32 * v + (1.0 - data.h) / 2.0 * (2.0 * PI as f32 * v).sin())
        .collect();

    let mut rho = vec![0.0f32; npts];
    let mut p = vec![0.0f32; npts];
    let mut bmag = vec![0.0f32; npts];
    let mut u = vec![FourVector::flat([0.0; 4]); npts];
    let mut b = vec![FourVector::flat([0.0; 4]); npts];

    for i in 0..npts {
        let zr = x0[i].data[1] as f32;
        let theta = x0[i].data[2] as f32;
        let (x1, x2) = transformbl2mksh(zr, theta, data.h);
        let mut lx1 =
            ((x1 - uniqx1[0]) / (uniqx1[nx1 - 1] - uniqx1[0]) * (nx1 as f32 - 1.0)).floor()
                as i32
                + 1;
        let mut lx2 =
            ((x2 - uniqx2[0]) / (uniqx2[nx2 - 1] - uniqx2[0]) * (nx2 as f32 - 1.0)).floor()
                as i32
                + 1;
        // upstream clamps lx2/ux2 only; lx1 is clamped here for memory
        // safety (out-of-grid rays read arbitrary in-bounds memory upstream)
        lx1 = lx1.clamp(1, nx1 as i32 - 1);
        lx2 = lx2.max(1);
        let ux1 = lx1 + 1;
        let mut ux2 = lx2 + 1;
        if ux2 > nx2 as i32 {
            ux2 = nx2 as i32;
        }
        // pole handling
        let dth = if ux2 != lx2 {
            uniqth[(ux2 - 1) as usize] - uniqth[(lx2 - 1) as usize]
        } else {
            uniqth[0]
        };
        let td = (theta - uniqth[(lx2 - 1) as usize]).abs() / dth;
        let mut rd = (zr - uniqr[(lx1 - 1) as usize])
            / (uniqr[(ux1 - 1) as usize] - uniqr[(lx1 - 1) as usize]);
        // inside the horizon-adjacent zone: nearest neighbour
        if uniqr[(lx1 - 1) as usize] <= 1.0 + (1.0 - a * a).sqrt() {
            rd = 1.0;
        }
        let x2l = lx2 - 1;
        let x2u = ux2 - 1;
        let x1l = (lx1 - 1) * nx2 as i32;
        let x1u = (ux1 - 1) * nx2 as i32;
        let idx = [
            (x1l + x2l) as usize,
            (x1l + x2u) as usize,
            (x1u + x2l) as usize,
            (x1u + x2u) as usize,
        ];
        let bil = |v: [f32; 4]| -> f32 {
            let w1 = v[0] * (1.0 - td) + v[1] * td;
            let w2 = v[2] * (1.0 - td) + v[3] * td;
            w1 * (1.0 - rd) + w2 * rd
        };
        let get4 = |arr: &[f32]| -> [f32; 4] {
            [arr[idx[0]], arr[idx[1]], arr[idx[2]], arr[idx[3]]]
        };
        let inside = x1 > uniqx1[0];
        let rho_v = bil(get4(&data.rho_arr));
        rho[i] = if inside { rho_v } else { 0.0 };
        let p_v = bil(get4(&data.p_arr));
        p[i] = if inside { p_v } else { 1.0 };
        let vrl0 = if inside { bil(get4(&data.vrl_arr)) } else { 0.0 };
        let vtl0 = if inside { bil(get4(&data.vtl_arr)) } else { 1.0 };
        let vpl0 = if inside { bil(get4(&data.vpl_arr)) } else { 0.0 };
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
        // metric for bmag (single-precision overload as gfortran resolves it)
        let metric = blmetric_cov_f32(zr, theta, a);
        b[i].assign_metric(metric);
        let bdotb = b[i].dot(&b[i]);
        bmag[i] = bdotb.max(0.0).sqrt() as f32;
        // four-velocity from the interpolated LNRF velocities
        let (vr0, vth0, vph0) = lnrf_frame_inv_f32(vrl0, vtl0, vpl0, zr, a, theta);
        u[i].data[1] = u[i].data[0] * vr0 as f64;
        u[i].data[2] = u[i].data[0] * vth0 as f64;
        u[i].data[3] = u[i].data[0] * vph0 as f64;
        let metric64 = blmetric_cov(zr as f64, theta as f64, a as f64);
        u[i].assign_metric(metric64);
    }
    (rho, p, bmag, u, b)
}

/// Upstream `get_harm_fluidvars`.
pub fn get_harm_fluidvars(data: &HarmData, x0: &[FourVector], a: f32, f: &mut crate::fluid::Fluid) {
    let (rho, p, bmag, u, b) = harm_vals(data, x0, a);
    f.rho = rho;
    f.p = p;
    f.bmag = bmag;
    f.u = u;
    f.b = b;
}

/// Upstream `scale_sim_units`: code units -> cgs.
#[allow(clippy::too_many_arguments)]
pub fn scale_sim_units(
    mbh: f64,
    mdotcgs: f64,
    mdot: f64,
    rho: &[f32],
    p: &[f32],
    bmag: &[f32],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = rho.len();
    let lcgs = GC * mbh * MSUN / C2;
    let tcgs = lcgs / grtrans_core::constants::C;
    let mut ncgs = vec![0.0f64; n];
    let mut bcgs = vec![0.0f64; n];
    let mut tempcgs = vec![0.0f64; n];
    for i in 0..n {
        let rhocgs = mdotcgs / mdot / lcgs.powi(3) * tcgs * (rho[i] as f64);
        ncgs[i] = rhocgs / MP;
        let pcgs = (p[i] as f64) * rhocgs / (rho[i] as f64) * C2;
        tempcgs[i] = pcgs / ncgs[i] / K;
        bcgs[i] = (bmag[i] as f64) * (rhocgs / (rho[i] as f64)).sqrt() * grtrans_core::constants::C;
        bcgs[i] *= (4.0 * PI).sqrt();
    }
    (ncgs, bcgs, tempcgs)
}
