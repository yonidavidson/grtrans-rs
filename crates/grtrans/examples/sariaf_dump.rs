//! Dump the SARIAF image for comparison.
use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use rayon::prelude::*;

fn main() {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/tmp/sariaf_ours.f64.bin".to_string());
    let args = initialize_pixels(
        true, 1, 0.6428, -0.5, 0.9375, 0.01, 1.0, 1.0, 2, -15.0, 15.0, -15.0, 15.0, 100, 100, 400,
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
    let mut loaded = load_fluid_model("SARIAF", 0.9375, &fargs);
    let opts = TraceOptions {
        ename: "POLSYNCHTH".to_string(),
        iname: "lsoda".to_string(),
        nvals: 1,
        freqs: vec![2.3e11],
        ..Default::default()
    };
    let sp = SourceParams {
        mbh: 6.2e8,
        mdot: 1.57e15,
        gminval: 1.0,
        muval: 0.25,
        ..Default::default()
    };
    let results: Vec<(usize, f64)> = (0..10000)
        .into_par_iter()
        .map(|i| {
            let (ray, status) = initialize_geodesic(&args, i);
            let mut ll = loaded.clone();
            let vals = trace_ray(&mut ll, &ray, status, &sp, &opts, false);
            (i, vals[0])
        })
        .collect();
    let _ = &mut loaded;
    let mut ours = vec![0.0f64; 10000];
    for (i, v) in results {
        ours[i] = v;
    }
    let mut bytes = Vec::new();
    for v in &ours {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(&out, bytes).unwrap();
    println!("wrote {out}");
}
