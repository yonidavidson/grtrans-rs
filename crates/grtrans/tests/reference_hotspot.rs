//! End-to-end reference test: HOTSPOT (Gaussian spot on a Keplerian orbit,
//! toroidal field) with POLSYNCHPL emissivity. Upstream tolerance: 1e-2.

use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;

#[test]
fn hotspot_matches_reference() {
    let spin = 0.9375;
    let mu0 = 0.6428;
    let nup = 400usize;
    let nvals = 4usize;
    let nu0 = 2.3e11;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -15.0, 15.0, -15.0, 15.0, 100, 100, nup,
    );
    let fargs = FluidArgs {
        rspot: 1.5,
        r0spot: 6.0,
        n0spot: 1e4,
        bl06: 1,
        mbh: 6.2e8,
        mdot: 1.57e15,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("HOTSPOT", spin, &fargs);
    let opts = TraceOptions {
        ename: "POLSYNCHPL".to_string(),
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
        p1: 3.5,
        p2: 3.5,
        gmax: 1e5,
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
    let full = format!(
        "{}/../../reference/fixtures/hotspot/ivals.f64.bin",
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
    println!("HOTSPOT relative error vs upstream reference: {err:e}");
    assert!(err < 1e-2, "HOTSPOT relative error {err:e} exceeds 1e-2");
}
