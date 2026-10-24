//! End-to-end reference test: SARIAF (semi-analytic RIAF).
//!
//! Fixture generated with the upstream Fortran binary (see
//! reference/fixtures/sariaf/inputs.in). Tolerance: upstream 1e-2.

use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;

fn fixture(path: &str) -> Vec<f64> {
    let full = format!(
        "{}/../../reference/fixtures/{path}",
        env!("CARGO_MANIFEST_DIR")
    );
    let bytes = std::fs::read(&full).unwrap_or_else(|_| panic!("missing fixture {full}"));
    bytes
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

#[test]
fn sariaf_matches_reference() {
    let spin = 0.9375;
    let mu0 = 0.6428;
    let nup = 400usize;
    let nvals = 1usize;
    let nu0 = 2.3e11;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -15.0, 15.0, -15.0, 15.0, 100, 100, nup,
    );
    let fargs = FluidArgs {
        nscl: 3e7,
        tscl: 1.7e11,
        nnthscl: 8e4,
        nnthp: 2.9,
        beta: 10.0,
        bl06: 0,
        mbh: 6.2e8,
        mdot: 1.57e15,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("SARIAF", spin, &fargs);
    let opts = TraceOptions {
        ename: "POLSYNCHTH".to_string(),
        iname: "lsoda".to_string(),
        nvals,
        freqs: vec![nu0],
        ..Default::default()
    };
    let sp = SourceParams {
        mbh: 6.2e8,
        mdot: 1.57e15,
        gminval: 1.0,
        muval: 0.25,
        ..Default::default()
    };
    let npix = 100 * 100;
    let results: Vec<(usize, Vec<f64>)> = (0..npix)
        .into_par_iter()
        .map(|i| {
            let (ray, status) = initialize_geodesic(&args, i);
            let mut ll = loaded.clone();
            let vals = trace_ray(&mut ll, &ray, status, &sp, &opts, false);
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
    let reference = fixture("sariaf/ivals.f64.bin");
    assert_eq!(reference.len(), ours.len());
    let num: f64 = ours
        .iter()
        .zip(reference.iter())
        .map(|(a, b)| (a - b).abs())
        .sum();
    let den: f64 = reference.iter().map(|b| b.abs()).sum();
    let err = num / den;
    println!("SARIAF relative error vs upstream reference: {err:e}");
    assert!(err < 1e-2, "SARIAF relative error {err:e} exceeds 1e-2");
}
