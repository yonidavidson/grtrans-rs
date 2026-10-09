//! Ray assembly: camera pixels to geodesic samples with BL coordinates,
//! wave vectors, affine parameters and turning-point parities.
//!
//! Direct translation of `geodesics.f90` (upstream GRTRANS):
//! `initialize_pixels`, `initialize_geodesic` and the per-pixel
//! `geokerr_wrapper`.

use crate::geokerr::camera::{initialize_camera_geokerr, CameraPixel};
use crate::geokerr::driver::geokerr;
use grtrans_core::four_vector::FourVector;
use grtrans_core::kerr::{blmetric_cov, calc_nullp};

/// Camera + geodesic arguments (upstream `geokerr_args`).
#[derive(Clone, Debug)]
pub struct GeokerrArgs {
    pub a: f64,
    pub mu0: f64,
    /// observer azimuth in units of pi (upstream multiplies by pi)
    pub phi0: f64,
    pub u0: f64,
    pub uout: f64,
    pub uin: f64,
    pub offset: f64,
    pub usegeor: i32,
    pub phit: i32,
    pub mufill: i32,
    pub kext: usize,
    pub next: usize,
    pub nup: usize,
    pub pixels: Vec<CameraPixel>,
    /// per-pixel initial coordinate time (zero unless nload > 1)
    pub t0: Vec<f64>,
}

/// One geodesic (upstream `geo`).
#[derive(Clone, Debug, Default)]
pub struct Ray {
    pub x: Vec<FourVector>,
    pub k: Vec<FourVector>,
    pub lambda: Vec<f64>,
    pub tpmarr: Vec<i32>,
    pub tprarr: Vec<i32>,
    pub npts: usize,
    // per-ray constants (upstream stores these in g%gk)
    pub alpha: f64,
    pub beta: f64,
    pub q2: f64,
    pub l: f64,
    pub su: f64,
    pub sm: f64,
    pub a: f64,
    pub mu0: f64,
    pub phi0: f64,
    pub t0: f64,
}

/// Upstream `initialize_pixels`: set up the camera and per-pixel data.
#[allow(clippy::too_many_arguments)]
pub fn initialize_pixels(
    use_geokerr: bool,
    standard: i32,
    mu0: f64,
    phi0: f64,
    a: f64,
    uout: f64,
    uin: f64,
    rcut: f64,
    nrotype: i32,
    a1: f64,
    a2: f64,
    b1: f64,
    b2: f64,
    nro: usize,
    nphi: usize,
    nup: usize,
) -> GeokerrArgs {
    let mut args = GeokerrArgs {
        a,
        mu0,
        phi0,
        u0: 0.0,
        uout,
        uin,
        offset: 0.0,
        usegeor: 0,
        phit: 1,
        mufill: 0,
        kext: 0,
        next: 0,
        nup,
        pixels: Vec::new(),
        t0: vec![0.0; nro * nphi],
    };
    if use_geokerr {
        if standard == 1 {
            args.usegeor = 0;
            args.phit = 1;
            if nup > 1 {
                args.mufill = 1;
                args.kext = (3.0 * nup as f64 / 400.0).ceil() as usize;
                args.next = (60.0 * nup as f64 / 400.0).ceil() as usize;
            } else {
                args.mufill = 0;
                args.kext = 0;
                args.next = 0;
            }
        } else {
            args.mufill = 0;
            args.phit = 1;
            args.usegeor = 1;
            args.kext = 0;
            args.next = 0;
        }
        let (u0, offset, pixels) = initialize_camera_geokerr(
            standard, a1, a2, b1, b2, rcut, nrotype, nro, nphi, nup, uout, mu0, a,
        );
        args.u0 = u0;
        args.offset = offset;
        args.pixels = pixels;
    }
    args
}

/// Upstream `initialize_geodesic`: run geokerr for pixel `i` and assemble
/// the ray. Returns `(ray, status)`; status = -1 marks an invalid ray
/// (negative radius sum), as upstream.
pub fn initialize_geodesic(args: &GeokerrArgs, i: usize) -> (Ray, i32) {
    let px = &args.pixels[i];
    let res = geokerr(
        args.u0,
        px.uf,
        args.uout,
        args.mu0,
        px.muf,
        args.a,
        px.l,
        px.q2,
        px.alpha,
        px.beta,
        px.tpm,
        px.tpr,
        px.su,
        px.sm,
        args.nup,
        args.offset,
        args.phit != 0,
        args.usegeor != 0,
        args.mufill != 0,
        args.kext,
        args.next,
    );

    let mut npts = res.nup;
    let mut i1 = 0usize;
    let mut i2 = npts;
    // Remove locations with invalid solutions (only for standard=2 rays).
    if args.usegeor == 1 && npts > 1 {
        if res.ufi[0] < 0.0 {
            let mut found = false;
            for ii in 1..npts {
                if res.ufi[ii] >= 0.0 {
                    i1 = ii;
                    found = true;
                    break;
                }
            }
            if !found {
                i1 = 0;
            }
        } else {
            i1 = 0;
        }
        if res.ufi[npts - 1] < 0.0 {
            for ii in 1..npts {
                if res.ufi[npts - ii - 1] >= 0.0 {
                    i2 = npts - ii;
                    break;
                }
            }
        } else {
            i2 = npts;
        }
        npts = i2 - i1;
    }

    let mut ray = Ray {
        npts,
        alpha: px.alpha,
        beta: px.beta,
        q2: px.q2,
        l: px.l,
        su: px.su,
        sm: px.sm,
        a: args.a,
        mu0: args.mu0,
        phi0: args.phi0,
        t0: args.t0[i],
        ..Default::default()
    };
    ray.x = Vec::with_capacity(npts);
    ray.k = Vec::with_capacity(npts);
    ray.lambda = vec![0.0; npts];
    ray.tpmarr = vec![0; npts];
    ray.tprarr = vec![0; npts];

    // upstream (geodesics.f90 line 218) computes acos(-1.) with a
    // single-precision literal, so phi0 is multiplied by the f32 value of
    // pi; reproduce that exactly.
    let phi0_pi = (std::f32::consts::PI) as f64 * args.phi0;
    for j in 0..npts {
        let idx = i1 + j;
        let r = 1.0 / res.ufi[idx];
        let theta = res.mufi[idx].acos();
        let mut phi = phi0_pi - res.dphi[idx];
        // pole-on viewing fix from upstream
        if args.mu0.abs() == 1.0 {
            let s = if args.mu0 >= 0.0 { 1.0 } else { -1.0 };
            phi += s * px.beta.atan2(px.alpha);
        }
        ray.x.push(FourVector::flat([0.0, r, theta, phi]));
    }
    if npts != 1 {
        for j in 0..npts {
            let idx = i1 + j;
            // upstream: LAMBDAI(g%npts) - LAMBDAI(i1:i2) with g%npts the
            // trimmed count (a Fortran index quirk reproduced exactly)
            ray.lambda[j] = res.lambdai[npts - 1] - res.lambdai[idx];
            ray.x[j].data[0] = res.dti[0] - res.dti[idx] - ray.t0;
            ray.tpmarr[j] = res.tpmi[idx];
            ray.tprarr[j] = res.tpri[idx];
        }
    } else {
        ray.lambda[0] = res.lambdai[0];
        ray.x[0].data[0] = res.dti[0];
        ray.tpmarr[0] = res.tpmi[0];
        ray.tprarr[0] = res.tpri[0];
    }

    // wave vector from the constants of motion
    for j in 0..npts {
        let r = ray.x[j].data[1];
        let mu = ray.x[j].data[2].cos();
        let su_signed = ray.su * if ray.tprarr[j] % 2 == 0 { 1.0 } else { -1.0 };
        let sm_signed = ray.sm * if ray.tpmarr[j] % 2 == 0 { 1.0 } else { -1.0 };
        let kv = calc_nullp(px.q2, px.l, args.a, r, mu, su_signed, sm_signed, true, true);
        ray.k.push(FourVector::flat(kv));
    }

    // assign BL covariant metrics (transpose of a symmetric matrix is self)
    for j in 0..npts {
        let metric = blmetric_cov(ray.x[j].data[1], ray.x[j].data[2], args.a);
        ray.x[j].assign_metric(metric);
        ray.k[j].assign_metric(metric);
    }

    let mut status = 1i32;
    let sum_r: f64 = ray.x.iter().map(|v| v.data[1]).sum();
    if sum_r < 0.0 {
        status = -1;
    }
    (ray, status)
}
