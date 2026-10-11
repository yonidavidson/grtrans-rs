//! Scan all POWERLAW rays serially with progress output to find slow pixels.
use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};
use std::io::Write;
use std::time::Instant;

fn main() {
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
    let iname = std::env::args().nth(2).unwrap_or_else(|| "lsoda".to_string());
    let opts = TraceOptions {
        ename: "POLSYNCHTH".to_string(), iname, nvals: 4,
        freqs: vec![3.45e11], ..Default::default()
    };
    let sp = SourceParams { mbh: 4e6, mdot: 1.5e15, gminval: 1.0, muval: 0.25, ..Default::default() };
    if std::env::args().nth(1).as_deref() == Some("dump") {
        let mut ours = vec![0.0f64; 40000 * 4];
        for i in 0..40000usize {
            let (ray, status) = initialize_geodesic(&args, i);
            let vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
            for q in 0..4 { ours[i * 4 + q] = vals[q]; }
        }
        let mut bytes = Vec::new();
        for v in &ours { bytes.extend_from_slice(&v.to_le_bytes()); }
        std::fs::write("/tmp/powerlaw_ours.f64.bin", bytes).unwrap();
        println!("wrote /tmp/powerlaw_ours.f64.bin");
        return;
    }
    let t0 = Instant::now();
    let mut worst = (0usize, 0.0f64);
    let lo: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(0);
    let hi: usize = std::env::args().nth(2).and_then(|v| v.parse().ok()).unwrap_or(40000);
    for i in lo..hi {
        println!("pixel {i}");
        std::io::stdout().flush().unwrap();
        let t1 = Instant::now();
        let (ray, status) = initialize_geodesic(&args, i);
        let dt_geo = t1.elapsed().as_secs_f64();
        println!("pixel {i}: geodesic done npts={} status={}", ray.npts, status);
        std::io::stdout().flush().unwrap();
        let t2 = Instant::now();
        let _vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
        let dt = t2.elapsed().as_secs_f64();
        if dt > worst.1 { worst = (i, dt); }
        if dt > 0.05 || dt_geo > 0.05 {
            println!("slow pixel {i}: geodesic {dt_geo:.3}s trace {dt:.3}s npts={}", ray.npts);
            std::io::stdout().flush().unwrap();
        }
        if i % 2000 == 0 {
            println!("progress {i} elapsed {:.1}s worst {} {:.3}s", t0.elapsed().as_secs_f64(), worst.0, worst.1);
            std::io::stdout().flush().unwrap();
        }
    }
    println!("done in {:.1}s, worst pixel {} {:.3}s", t0.elapsed().as_secs_f64(), worst.0, worst.1);
}
