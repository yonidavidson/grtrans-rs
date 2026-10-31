//! End-to-end reference test: POWERLAW fluid with MAXJUTT emissivity
//! (6-weight Maxwell-Jüttner decomposition, delta=2). Exercises the
//! `otherargs` plumbing through the driver. Upstream tolerance: 2e-2.

use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;

#[test]
fn powerlaw_maxjutt_matches_reference() {
    let spin = 0.0;
    let mu0 = 0.5;
    let nx = 200usize;
    let ny = 200usize;
    let nup = 1600usize;
    let nvals = 4usize;
    let nu0 = 3.45e11;
    let args = initialize_pixels(
        true,
        1,
        mu0,
        -0.5,
        spin,
        5e-5,
        1.0,
        4000.0,
        1,
        1200.0,
        4000.0,
        0.0,
        2.0 * std::f64::consts::PI,
        nx,
        ny,
        nup,
    );
    let fargs = FluidArgs {
        nscl: 1e5,
        tscl: 5e11,
        nnthscl: 8e4,
        nnthp: 2.9,
        beta: 10.0,
        np: 0.0,
        tp: 0.0,
        rin: 3200.0,
        rout: 3300.0,
        thin: -0.02,
        thout: 0.02,
        phiin: 0.0,
        phiout: 1e4,
        mbh: 4e6,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("POWERLAW", spin, &fargs);
    let opts = TraceOptions {
        ename: "MAXJUTT".to_string(),
        iname: "lsoda".to_string(),
        nvals,
        freqs: vec![nu0],
        otherargs: vec![
            2.0, 0.10714295, 0.4459917, 1.34219997, 2.17450422, 1.36515333, 0.52892829,
        ],
        ..Default::default()
    };
    let sp = SourceParams {
        mbh: 4e6,
        mdot: 1.5e15,
        gminval: 1.0,
        muval: 0.25,
        ..Default::default()
    };
    let npix = nx * ny;
    let results: Vec<(usize, Vec<f64>)> = (0..npix)
        .into_par_iter()
        .map(|i| {
            let (ray, status) = initialize_geodesic(&args, i);
            let mut loaded_local = loaded.clone();
            let vals = trace_ray(&mut loaded_local, &ray, status, &sp, &opts, false);
            (i, vals)
        })
        .collect();
    let _ = &mut loaded;
    let mut ours = vec![0.0f64; npix * nvals];
    for (i, vals) in results {
        for q in 0..nvals {
            ours[i * nvals + q] = vals[q];
        }
    }
    let full = format!(
        "{}/../../reference/fixtures/powerlaw_maxjutt/ivals.f64.bin",
        env!("CARGO_MANIFEST_DIR")
    );
    let bytes = std::fs::read(&full).unwrap();
    let reference: Vec<f64> = bytes
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    assert_eq!(reference.len(), ours.len());
    let num: f64 = ours
        .iter()
        .zip(reference.iter())
        .map(|(a, b)| (a - b).abs())
        .sum();
    let den: f64 = reference.iter().map(|b| b.abs()).sum();
    let err = num / den;
    println!("POWERLAW MAXJUTT relative error vs upstream reference: {err:e}");
    assert!(err < 2e-2, "MAXJUTT relative error {err:e} exceeds 2e-2");
}
