//! Dump the THINDISK image to a raw file for comparison tooling.
//! Run: cargo run --release -p grtrans --example thindisk_dump -- /tmp/ours.f64.bin

use grtrans::driver::{trace_ray, SourceParams, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs};

fn main() {
    let out_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/tmp/thindisk_ours.f64.bin".to_string());
    let spin = 0.9;
    let mu0 = 0.26;
    let nx = 100usize;
    let ny = 100usize;
    let nfreq = 25usize;
    let nvals = 4usize;
    let fmin: f64 = 2.41e16;
    let fmax: f64 = 6.31e18;
    let args = initialize_pixels(
        true, 2, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -21.0, 21.0, -21.0, 21.0, nx, ny, 1,
    );
    let fargs = FluidArgs {
        mdot: 0.1,
        mbh: 10.0,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("THINDISK", spin, &fargs);
    let freqs: Vec<f64> = (0..nfreq)
        .map(|i| fmin * ((i as f64) * (fmax / fmin).ln() / ((nfreq - 1) as f64)).exp())
        .collect();
    let opts = TraceOptions {
        ename: "BBPOL".to_string(),
        iname: "lsoda".to_string(),
        nvals,
        freqs,
    };
    let sp = SourceParams {
        mdot: 1.5e15,
        mbh: 10.0,
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
    let mut bytes = Vec::with_capacity(ours.len() * 8);
    for v in &ours {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(&out_path, bytes).unwrap();
    println!("wrote {out_path}");
}
