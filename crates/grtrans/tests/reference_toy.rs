//! End-to-end self-consistency test: TOY fluid model (code-comparison toy)
//! with POLSYNCHPL emissivity and the documented aleph/bb choices (upstream's
//! aleph uses stale f%u and bb is unassigned, so the upstream image is
//! garbage and not a valid reference). toy_vals and f%u are validated in
//! toy_fixture.rs; here we check the full path is finite and physical.

use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;

#[test]
fn toy_interp_self_consistent() {
    let spin = 0.9375;
    let mu0 = 0.6428;
    let nup = 400usize;
    let nvals = 4usize;
    let nu0 = 2.3e11;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -15.0, 15.0, -15.0, 15.0, 100, 100, nup,
    );
    let fargs = FluidArgs {
        nscl: 3e7,
        tscl: 1.7,
        beta: 2.0,
        bl06: 0,
        mbh: 6.2e8,
        mdot: 1.57e15,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("TOY", spin, &fargs);
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
    let nfinite = ours.iter().filter(|v| v.is_finite()).count();
    assert_eq!(nfinite, ours.len(), "non-finite intensities present");
    let nonzero = ours.iter().filter(|v| **v != 0.0).count();
    let max = ours.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!(
        "TOY self-consistency: nonzero {nonzero}/{} max {max:e}",
        ours.len()
    );
    assert!(nonzero > 0, "toy should emit somewhere");
    assert!(max.is_finite(), "unphysical maximum");
}
