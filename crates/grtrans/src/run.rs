//! Image orchestration: camera setup, parallel pixel loop, binary output.
//!
//! Mirrors `pgrtrans.f90` / `grtrans.f90` for the single-time-step case
//! (the time-dependent fluid advance is not ported).

use crate::driver::{trace_ray, TraceOptions};
use crate::inputs::{freqs, mdots, mus, Inputs};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels, GeokerrArgs};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

/// One camera image (all pixels for one mu/mdot/frequency combination).
#[derive(Clone, Debug, Default)]
pub struct Image {
    pub nx: usize,
    pub ny: usize,
    pub nvals: usize,
    pub nu: f64,
    /// pixel locations (alpha, beta), one per pixel
    pub ab: Vec<[f32; 2]>,
    /// intensities, row-major [pixel][stokes]
    pub ivals: Vec<f32>,
}

/// Run the full pipeline described by `inputs`.
pub fn run(inp: &Inputs) -> Vec<Image> {
    let nro = inp.nn[0] as usize;
    let nphi = inp.nn[1] as usize;
    let nup = inp.nn[2] as usize;
    let npix = nro * nphi;
    let nvals = inp.nvals as usize;

    // i1/i2 handling (upstream grtrans_batch clamps them)
    let mut i1 = if inp.i1 < 0 { 1 } else { inp.i1 };
    let mut i2 = if inp.i2 < 0 { npix as i64 } else { inp.i2 };
    if i1 > npix as i64 - 1 {
        i1 = npix as i64 - 1;
    }
    if i2 > npix as i64 {
        i2 = npix as i64;
    }

    let fargs = FluidArgs {
        dfile: inp.fdfile.clone(),
        hfile: inp.fhfile.clone(),
        gfile: inp.fgfile.clone(),
        sim: inp.fsim.clone(),
        nt: inp.fnt as i32,
        indf: inp.findf as i32,
        nfiles: inp.fnfiles as i32,
        jonfix: inp.fjonfix as i32,
        offset: inp.foffset as i32,
        dindf: inp.fdindf as i32,
        magcrit: inp.fmagcrit as i32,
        scalefac: inp.fscalefac,
        mdot: inp.fmdot,
        mbh: inp.mbh,
        rin: inp.frin,
        rout: inp.frout,
        sigcut: inp.sigcut,
        tscl: inp.ftscl,
        rscl: inp.frscl,
        nw: inp.fnw as i32,
        nfreq_tab: inp.fnfreq_tab as i32,
        nr: inp.fnr as i32,
        wmin: inp.fwmin,
        wmax: inp.fwmax,
        fmin: inp.ffmin,
        fmax: inp.ffmax,
        rmax: inp.frmax,
        sigt: inp.fsigt,
        fcol: inp.ffcol,
        rspot: inp.frspot,
        r0spot: inp.fr0spot,
        n0spot: inp.fn0spot,
        nscl: inp.fnscl,
        nnthscl: inp.fnnthscl,
        nnthp: inp.fnnthp,
        beta: inp.fbeta,
        bl06: inp.fbl06 as i32,
        np: inp.fnp,
        tp: inp.ftp,
        thin: inp.fthin,
        thout: inp.fthout,
        phiin: inp.fphiin,
        phiout: inp.fphiout,
        ..Default::default()
    };
    let mut loaded = load_fluid_model(&inp.fname, inp.spin, &fargs);

    let freq_grid = freqs(inp);
    let mdot_grid = mdots(inp);
    let mu_grid = mus(inp);
    let opts = TraceOptions {
        ename: inp.ename.clone(),
        iname: inp.iname.clone(),
        nvals,
        freqs: freq_grid.clone(),
        ..Default::default()
    };

    let mut images = Vec::new();
    for &mu0 in &mu_grid {
        let args: GeokerrArgs = initialize_pixels(
            inp.use_geokerr,
            inp.standard,
            mu0,
            inp.phi0,
            inp.spin,
            inp.uout,
            inp.uin,
            inp.rcut,
            inp.nrotype,
            inp.gridvals[0],
            inp.gridvals[1],
            inp.gridvals[2],
            inp.gridvals[3],
            nro,
            nphi,
            nup,
        );
        for &mdot in &mdot_grid {
            let sp = SourceParams {
                nfac: inp.ftscl,
                bfac: inp.frscl,
                mbh: inp.mbh,
                mdot,
                p1: inp.p1,
                p2: inp.p2,
                gmax: inp.gmax,
                gminval: inp.gmin,
                jetalphaval: inp.jetalpha,
                muval: inp.muval,
                sigcut: inp.sigcut,
                gmin: vec![inp.gmin; nup],
                mu: vec![inp.muval; nup],
                jetalpha: vec![inp.jetalpha; nup],
                type_: 0,
            };
            // parallel over pixels (each ray is independent; results are
            // deterministic and order-preserving)
            let results: Vec<(usize, Vec<f64>)> = (i1..=i2)
                .into_par_iter()
                .map(|i| {
                    let (ray, status) = initialize_geodesic(&args, (i - 1) as usize);
                    let mut loaded_local = loaded.clone();
                    let vals = trace_ray(&mut loaded_local, &ray, status, &sp, &opts, false);
                    ((i - 1) as usize, vals)
                })
                .collect();
            // merge the per-model state mutations (THINDISK's rin update is
            // the only one and is idempotent across identical rays)
            let _ = &mut loaded;

            for k in 0..freq_grid.len() {
                let mut img = Image {
                    nx: nro,
                    ny: nphi,
                    nvals,
                    nu: freq_grid[k],
                    ab: vec![[0.0; 2]; npix],
                    ivals: vec![0.0f32; npix * nvals],
                };
                for (i, p) in args.pixels.iter().enumerate() {
                    img.ab[i] = [p.alpha as f32, p.beta as f32];
                }
                for (i, vals) in &results {
                    for q in 0..nvals {
                        img.ivals[i * nvals + q] = vals[k * nvals + q] as f32;
                    }
                }
                images.push(img);
            }
        }
    }
    images
}

/// Write one image in the upstream plain-binary format
/// (`kwrite_raytrace_camera`, cflag=0): Fortran unformatted records
/// `nx,ny,nvals`; `nkey`; one record per key; `ab`; `ivals`.
pub fn write_binary(images: &[Image], path: &Path) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    for (idx, img) in images.iter().enumerate() {
        let mut w = |bytes: &[u8]| -> std::io::Result<()> {
            f.write_all(&(bytes.len() as i32).to_le_bytes())?;
            f.write_all(bytes)?;
            f.write_all(&(bytes.len() as i32).to_le_bytes())?;
            Ok(())
        };
        let mut rec = Vec::new();
        for v in [img.nx as i32, img.ny as i32, img.nvals as i32] {
            rec.extend_from_slice(&v.to_le_bytes());
        }
        w(&rec)?;
        // nkey: upstream writes 3 keys (nx, ny, nu) via the FITS path only;
        // the binary path writes nkey=0 when not supplied
        let _ = idx;
        w(&0i32.to_le_bytes())?;
        let mut ab = Vec::with_capacity(img.ab.len() * 8);
        for p in &img.ab {
            ab.extend_from_slice(&p[0].to_le_bytes());
            ab.extend_from_slice(&p[1].to_le_bytes());
        }
        w(&ab)?;
        let mut iv = Vec::with_capacity(img.ivals.len() * 4);
        for v in &img.ivals {
            iv.extend_from_slice(&v.to_le_bytes());
        }
        w(&iv)?;
    }
    Ok(())
}
