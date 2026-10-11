//! Timing benchmark: trace N POWERLAW rays.
use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use std::time::Instant;

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(8);
    let args = initialize_pixels(
        true, 1, 0.5, -0.5, 0.0, 5e-5, 1.0, 4000.0, 1, 1200.0, 4000.0, 0.0,
        2.0 * std::f64::consts::PI, 200, 200, 1600,
    );
    let fargs = FluidArgs {
        nscl: 1e5, tscl: 5e11, nnthscl: 8e4, nnthp: 2.9, beta: 10.0,
        rin: 3200.0, rout: 3300.0, thin: -0.02, thout: 0.02, mbh: 4e6,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("POWERLAW", 0.0, &fargs);
    let opts = TraceOptions {
        ename: "POLSYNCHTH".to_string(), iname: "lsoda".to_string(), nvals: 4,
        freqs: vec![3.45e11], ..Default::default()
    };
    let sp = SourceParams { mbh: 4e6, mdot: 1.5e15, gminval: 1.0, muval: 0.25, ..Default::default() };
    let t0 = Instant::now();
    let stride = 40000 / n.max(1);
    let mut worst = (0usize, 0.0f64);
    let mut total = 0.0f64;
    let mut nz = 0usize;
    for k in 0..n {
        let i = k * stride;
        let (ray, status) = initialize_geodesic(&args, i);
        let t1 = Instant::now();
        let vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
        let dt = t1.elapsed().as_secs_f64();
        total += dt;
        if dt > worst.1 { worst = (i, dt); }
        if vals.iter().any(|v| *v != 0.0) { nz += 1; }
        if k < 12 || dt > 0.5 {
            println!("pixel {i}: npts={} trace {dt:.3}s I0={:e}", ray.npts, vals[0]);
        }
    }
    println!("total {total:.3}s for {n} rays; worst pixel {} {:.3}s; nonzero {nz}", worst.0, worst.1);
}
