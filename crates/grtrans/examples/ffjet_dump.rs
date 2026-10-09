//! Dump the FFJET image to a raw file for comparison tooling.
use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{load_fluid_model, FluidArgs, SourceParams};

fn main() {
    let out_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/tmp/ffjet_ours.f64.bin".to_string());
    let spin = 0.998;
    let mu0 = 0.906;
    let nx = 100usize;
    let ny = 100usize;
    let nup = 400usize;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -40.0, 20.0, -20.0, 40.0, nx, ny, nup,
    );
    let dump = format!(
        "{}/../../reference/fixtures/ffjet/m87bl09rfp10xi5a998fluidvars.bin",
        env!("CARGO_MANIFEST_DIR")
    );
    let fargs = FluidArgs {
        dfile: dump,
        mdot: 1.5e15,
        mbh: 3.4e9,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("FFJET", spin, &fargs);
    let mut sp = SourceParams {
        nfac: 2.0,
        bfac: 70.0,
        mbh: 3.4e9,
        mdot: 1.5e15,
        p1: 3.5,
        p2: 3.5,
        gmax: 1e5,
        gminval: 100.0,
        jetalphaval: 0.02,
        muval: 0.25,
        sigcut: 1e10,
        ..Default::default()
    };
    sp.gmin = vec![100.0; nup];
    sp.mu = vec![0.25; nup];
    sp.jetalpha = vec![0.02; nup];
    let opts = TraceOptions {
        ename: "POLSYNCHPL".to_string(),
        iname: std::env::args()
            .nth(2)
            .unwrap_or_else(|| "lsoda".to_string()),
        nvals: 4,
        freqs: vec![3.45e11],
        ..Default::default()
    };
    let npix = nx * ny;
    let mut ours = vec![0.0f64; npix * 4];
    for i in 0..npix {
        let (ray, status) = initialize_geodesic(&args, i);
        let vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
        for q in 0..4 {
            ours[i * 4 + q] = vals[q];
        }
    }
    let mut bytes = Vec::with_capacity(ours.len() * 8);
    for v in &ours {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(&out_path, bytes).unwrap();
    println!("wrote {out_path}");
}
