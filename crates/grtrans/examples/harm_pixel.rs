//! Print the HARM chain summary for one pixel (0-based index from argv).
use grtrans::driver::TraceOptions;
use grtrans_core::kerr::comoving_ortho;
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::emissivity::{invariant_emis_g2, Emis};
use grtrans_physics::fluid::{convert_fluid_vars, get_fluid_vars, initialize_fluid_model, load_fluid_model, FluidArgs, SourceParams};

fn main() {
    let pix: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(5774);
    let args = initialize_pixels(
        true, 1, 0.6428, -0.5, 0.9375, 0.04, 1.0, 1.0, 2, -13.0, 13.0, -13.0, 13.0, 150, 150, 400,
    );
    let dump = format!("{}/../../reference/fixtures/harm/dump", env!("CARGO_MANIFEST_DIR"));
    let fargs = FluidArgs {
        dfile: dump, hfile: "dump040".to_string(), indf: 40, nt: 1,
        mbh: 4e6, mdot: 1.57e15, ..Default::default()
    };
    let mut loaded = load_fluid_model("HARM", 0.9375, &fargs);
    let sp = SourceParams { mbh: 4e6, mdot: 1.57e15, gminval: 1.0, muval: 0.25, ..Default::default() };
    let (ray, status) = initialize_geodesic(&args, pix);
    let n = ray.npts;
    println!("pixel {pix}: status={status} npts={n} alpha={:.4} beta={:.4}", ray.alpha, ray.beta);
    let mut f = initialize_fluid_model(&loaded, n);
    get_fluid_vars(&mut loaded, &ray.x, &ray.k, 0.9375, &mut f);
    let (ncgs, _ncgsnth, bcgs, tcgs) = convert_fluid_vars(&f, &sp);
    let imax = (0..n).max_by(|a, b| ncgs[*a].partial_cmp(&ncgs[*b]).unwrap()).unwrap();
    println!("  ncgs max {:.6e} at {imax} | tcgs max {:.6e} | bcgs max {:.6e}",
        ncgs[imax], tcgs.iter().cloned().fold(f64::MIN, f64::max), bcgs.iter().cloned().fold(f64::MIN, f64::max));
    // emissivity chain
    let mut u = f.u.clone(); let mut b = f.b.clone(); let mut kv = ray.k.clone();
    let mut rshift = vec![0.0f64; n]; let mut ang = vec![0.0f64; n]; let mut cosne = vec![0.0f64; n];
    let mut s2xi = vec![0.0f64; n]; let mut c2xi = vec![0.0f64; n];
    for i in 0..n {
        let co = comoving_ortho(ray.x[i].data[1], ray.x[i].data[2], 0.9375, ray.alpha, ray.beta, ray.mu0, &mut u[i], &mut b[i], &mut kv[i]);
        rshift[i] = co.g; ang[i] = co.ang; cosne[i] = co.cosne; s2xi[i] = co.s2xi; c2xi[i] = co.c2xi;
    }
    let mut e = Emis::select("POLSYNCHTH");
    e.initialize(n, &rshift, &ang, &cosne);
    e.assign_synch_params(&ncgs, &_ncgsnth, &bcgs, &tcgs);
    let nu: Vec<f64> = rshift.iter().map(|g| 2.3e11 / g).collect();
    e.calc_emissivity(&nu);
    let fac: f64 = (0..n).map(|i| e.j[i*4] * (ray.lambda[0] - ray.lambda[i])).sum();
    let lbh = grtrans_core::constants::G * 4e6 * grtrans_core::constants::MSUN / grtrans_core::constants::C2;
    let mut j1 = vec![0.0f64; n]; let mut k1 = vec![0.0f64; n];
    for i in 0..n {
        j1[i] = e.j[i*4] * rshift[i] * rshift[i] / fac;
        k1[i] = e.kcoef[i*e.nk] * lbh / rshift[i];
    }
    let jm = j1.iter().cloned().fold(f64::MIN, f64::max);
    let jmi = j1.iter().enumerate().max_by(|a,b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
    let km = k1.iter().cloned().fold(f64::MIN, f64::max);
    let kmi = k1.iter().enumerate().max_by(|a,b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
    println!("  g max {:.6e} | j1 max {jm:.6e} at {jmi} | K1 max {km:.6e} at {kmi} | fac {fac:.6e}", rshift.iter().cloned().fold(f64::MIN, f64::max));
    let opts = TraceOptions { ename: "POLSYNCHTH".to_string(), iname: "lsoda".to_string(), nvals: 1, freqs: vec![2.3e11], ..Default::default() };
    let vals = grtrans::driver::trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
    println!("  I = {:.6e}", vals[0]);
    // print tau max via the same path
    let tau: Vec<f64> = grtrans_transfer::calc_opt_depth(&ray.lambda, &k1).iter().map(|v| -v).collect();
    println!("  tau max {:.6e}", tau.iter().cloned().fold(f64::MIN, f64::max));
    // replicate the driver integration and print the profile
    let mut e2 = Emis::select("POLSYNCHTH");
    e2.initialize(n, &rshift, &ang, &cosne);
    e2.assign_synch_params(&ncgs, &_ncgsnth, &bcgs, &tcgs);
    e2.calc_emissivity(&nu);
    let fac2: f64 = (0..n).map(|i| e2.j[i*4] * (ray.lambda[0] - ray.lambda[i])).sum();
    grtrans_physics::emissivity::rotate_emis(&mut e2, &s2xi, &c2xi);
    let b194 = 194*e2.nk;
    println!("   post-rot idx 194: j={:?}", &e2.j[194*4..194*4+4]);
    println!("   post-rot idx 194: K={:?}", &e2.kcoef[b194..b194+7]);
    invariant_emis_g2(&mut e2, &rshift);
    for v in e2.j.iter_mut() { *v /= fac2; }
    for v in e2.kcoef.iter_mut() { *v *= lbh; }
    let jv: Vec<[f64;4]> = (0..n).map(|i| [e2.j[i*4], e2.j[i*4+1], e2.j[i*4+2], e2.j[i*4+3]]).collect();
    let kv: Vec<[f64;7]> = (0..n).map(|i| { let b=i*e2.nk; [e2.kcoef[b],e2.kcoef[b+1],e2.kcoef[b+2],e2.kcoef[b+3],e2.kcoef[b+4],e2.kcoef[b+5],e2.kcoef[b+6]] }).collect();
    let tau2: Vec<f64> = grtrans_transfer::calc_opt_depth(&ray.lambda, &kv.iter().map(|k| k[0]).collect::<Vec<_>>()).iter().map(|v| -v).collect();
    {
        let mut bytes = Vec::new();
        for i in 0..n {
            for c in 0..4 { bytes.extend_from_slice(&jv[i][c].to_le_bytes()); }
            for c in 0..7 { bytes.extend_from_slice(&kv[i][c].to_le_bytes()); }
            bytes.extend_from_slice(&ray.lambda[i].to_le_bytes());
        }
        std::fs::write("/tmp/harm_arrays.f64.bin", bytes).unwrap();
    }
    let (inten, nptsout) = grtrans_transfer::integrate(grtrans_transfer::Method::Lsoda, &ray.lambda, &jv, &kv, &tau2, 1e-2, 0.1, 1e-8, 1e-6);
    let (inten_tight, _) = grtrans_transfer::integrate(grtrans_transfer::Method::Lsoda, &ray.lambda, &jv, &kv, &tau2, 1e-2, 0.01, 1e-12, 1e-10);
    let (inten_tight2, _) = grtrans_transfer::integrate(grtrans_transfer::Method::Lsoda, &ray.lambda, &jv, &kv, &tau2, 1e-2, 0.001, 1e-14, 1e-12);
    let (inten_delo, _) = grtrans_transfer::integrate(grtrans_transfer::Method::Delo, &ray.lambda, &jv, &kv, &tau2, 1e-2, 0.1, 1e-8, 1e-6);
    let (inten_formal, _) = grtrans_transfer::integrate(grtrans_transfer::Method::Formal, &ray.lambda, &jv, &kv, &tau2, 1e-2, 0.1, 1e-8, 1e-6);
    println!("  nptsout={nptsout}");
    for i in [194usize] {
        let b=i*e2.nk;
        println!("   full idx {i}: j={:?}", &e2.j[i*4..i*4+4]);
        println!("   full idx {i}: K={:?}", &e2.kcoef[b..b+7]);
        println!("   full idx {i}: s2xi={:e} c2xi={:e} ang={:e} cosne={:e}", s2xi[i], c2xi[i], ang[i], cosne[i]);
    }
    for i in [170usize, 180, 190, 194, 200, 210, 250, 260, 270] {
        println!("   delo-profile idx {i}: delo={:e} formal={:e} lsoda={:e}", inten_delo[i][0], inten_formal[i][0], inten[i][0]);
    }
    for i in [0usize, 50, 100, 150, 194, 250, 350, 400, 453] {
        println!("   idx {i}: lam={:e} I={:e} tight={:e} tighter={:e} delo={:e} formal={:e}", ray.lambda[i], inten[i][0], inten_tight[i][0], inten_tight2[i][0], inten_delo[i][0], inten_formal[i][0]);
    }
}
