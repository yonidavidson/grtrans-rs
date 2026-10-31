//! End-to-end test: PHATDISK (Novikov-Thorne multicolor thin disk) with the
//! INTERP emissivity, as a self-consistency check.
//!
//! Why not an upstream-image comparison: upstream `interpemis_noabs` only
//! sets the emission column and leaves the absorption columns of the K buffer
//! UNINITIALIZED. The driver then uses that garbage as the absorption. The
//! geodebug for one pixel shows K1 ~ 1e5-1e6 and tau ~ 1.5e8 (optically thick
//! from garbage memory), so the upstream PHATDISK image depends on
//! uninitialized memory and is not a valid reference. We implement the
//! intended no-absorption behavior (K=0), validated at the component level
//! (test_phatdisk_fluid and test_interp fixtures).

use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;

#[test]
fn phatdisk_interp_self_consistent() {
    let spin = 0.9375;
    let mu0 = 0.6428;
    let nup = 400usize;
    let nvals = 1usize;
    let nu0 = 2.3e11;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -15.0, 15.0, -15.0, 15.0, 100, 100, nup,
    );
    let fargs = FluidArgs {
        nw: 500,
        wmin: 1e-4,
        wmax: 1e4,
        nfreq_tab: 100,
        fmin: 1e10,
        fmax: 1e15,
        rmax: 1e4,
        nr: 500,
        sigt: 0.4,
        fcol: 1.7,
        rin: -1.0,
        rout: 1e8,
        mbh: 6.2e8,
        mdot: 1.57e15,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("PHATDISK", spin, &fargs);
    let opts = TraceOptions {
        ename: "INTERP".to_string(),
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
    // No-absorption model: intensity must be finite and non-negative, and the
    // disk (which most rays cross) should emit at most pixels.
    let nfinite = ours.iter().filter(|v| v.is_finite()).count();
    assert_eq!(nfinite, ours.len(), "non-finite intensities present");
    let nneg = ours.iter().filter(|v| **v < 0.0).count();
    assert_eq!(nneg, 0, "negative intensity in a no-absorption model");
    let nonzero = ours.iter().filter(|v| **v != 0.0).count();
    let max = ours.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!(
        "PHATDISK self-consistency: nonzero {nonzero}/{} max {max:e}",
        ours.len()
    );
    assert!(nonzero > npix / 2, "disk should emit at most pixels");
    assert!(max > 0.0 && max.is_finite(), "unphysical maximum");
}
