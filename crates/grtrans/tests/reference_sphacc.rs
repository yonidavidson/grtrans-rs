//! End-to-end reference test: SPHACC (Bondi/Michel spherical accretion).
//!
//! Reproduces the upstream regression case
//! (run_grtrans_test_problems_public.py) from reference/fixtures/sphacc.
//! Upstream tolerance: 1e-1 relative (profile and spectrum).

use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};

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
fn sphacc_matches_reference() {
    let spin = 0.0;
    let mu0 = 0.1;
    let nx = 10000usize;
    let ny = 1usize;
    let nup = 100usize;
    let nfreq = 25usize;
    let nvals = 1usize;
    let fmin: f64 = 1e8;
    let fmax: f64 = 1e15;

    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.0025, 1.0, 1.0, 2, 0.0, 400.0, 0.0, 0.0, nx, ny, nup,
    );
    let fargs = FluidArgs {
        mbh: 1.0,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("SPHACC", spin, &fargs);

    let freqs: Vec<f64> = (0..nfreq)
        .map(|i| fmin * (((i as f64) * (fmax / fmin).ln()) / ((nfreq - 1) as f64)).exp())
        .collect();
    let opts = TraceOptions {
        ename: "SYNCHTHAV".to_string(),
        iname: "lsoda".to_string(),
        nvals,
        freqs,
        ..Default::default()
    };
    let sp = SourceParams {
        mbh: 1.0,
        mdot: 1.0,
        ..Default::default()
    };

    let npix = nx * ny;
    let mut ours = vec![0.0f64; npix * nvals * nfreq];
    for i in 0..npix {
        let (ray, status) = initialize_geodesic(&args, i);
        let vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
        for k in 0..nfreq {
            for q in 0..nvals {
                ours[(i * nvals + q) * nfreq + k] = vals[k * nvals + q];
            }
        }
    }
    // upstream spectrum for ny==1 (grtrans_batch.calc_spec):
    // spec = sum_pixels I * alpha * da * 2*pi with da = ab[1,0]-ab[0,0]
    let da = args.pixels[1].alpha - args.pixels[0].alpha;
    let mut spectrum = vec![0.0f64; nfreq];
    for k in 0..nfreq {
        spectrum[k] = (0..npix)
            .map(|i| ours[i * nvals * nfreq + k] * args.pixels[i].alpha)
            .sum::<f64>()
            * da
            * 2.0
            * std::f64::consts::PI;
    }

    let reference = fixture("sphacc/ivals.f64.bin");
    assert_eq!(reference.len(), ours.len());
    let num: f64 = ours
        .iter()
        .zip(reference.iter())
        .map(|(a, b)| (a - b).abs())
        .sum();
    let den: f64 = reference.iter().map(|b| b.abs()).sum();
    let err = num / den;
    println!("SPHACC image relative error vs upstream reference: {err:e}");
    assert!(
        err < 1e-1,
        "SPHACC image relative error {err:e} exceeds 1e-1"
    );

    let spec_ref = fixture("sphacc/spec.f64.bin");
    let spec_err = spectrum
        .iter()
        .zip(spec_ref.iter())
        .map(|(a, b)| (a - b).abs())
        .sum::<f64>()
        / spec_ref.iter().map(|b| b.abs()).sum::<f64>();
    println!("SPHACC spectrum relative error vs upstream reference: {spec_err:e}");
    assert!(
        spec_err < 1e-1,
        "SPHACC spectrum error {spec_err:e} exceeds 1e-1"
    );
}
