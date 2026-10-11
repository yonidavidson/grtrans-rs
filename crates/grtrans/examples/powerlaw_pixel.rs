//! Print the POWERLAW chain summary for one pixel (0-based index from argv).
use grtrans::driver::TraceOptions;
use grtrans_core::kerr::comoving_ortho;
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::emissivity::{invariant_emis_g2, Emis};
use grtrans_physics::fluid::{convert_fluid_vars, get_fluid_vars, initialize_fluid_model, load_fluid_model, FluidArgs, SourceParams};

fn main() {
    let pix: usize = std::env::args().nth(1).and_then(|v| v.parse().ok()).unwrap_or(10042);
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
    let sp = SourceParams { mbh: 4e6, mdot: 1.5e15, gminval: 1.0, muval: 0.25, ..Default::default() };
    let (ray, status) = initialize_geodesic(&args, pix);
    let n = ray.npts;
    println!("pixel {pix}: status={status} npts={n} alpha={:.4} beta={:.4}", ray.alpha, ray.beta);
    let mut f = initialize_fluid_model(&loaded, n);
    get_fluid_vars(&mut loaded, &ray.x, &ray.k, 0.0, &mut f);
    let (ncgs, _ncgsnth, bcgs, tcgs) = convert_fluid_vars(&f, &sp);
    let nmax = ncgs.iter().cloned().fold(f64::MIN, f64::max);
    let tmax = tcgs.iter().cloned().fold(f64::MIN, f64::max);
    let bmax = bcgs.iter().cloned().fold(f64::MIN, f64::max);
    println!("  n max {nmax:e} T max {tmax:e} B max {bmax:e}");
    // emissivity chain
    let mut u = f.u.clone(); let mut b = f.b.clone(); let mut kv = ray.k.clone();
    let mut rshift = vec![0.0f64; n]; let mut ang = vec![0.0f64; n]; let mut cosne = vec![0.0f64; n];
    for i in 0..n {
        let co = comoving_ortho(ray.x[i].data[1], ray.x[i].data[2], 0.0, ray.alpha, ray.beta, ray.mu0, &mut u[i], &mut b[i], &mut kv[i]);
        rshift[i] = co.g; ang[i] = co.ang; cosne[i] = co.cosne;
    }
    let mut e = Emis::select("POLSYNCHTH");
    e.initialize(n, &rshift, &ang, &cosne);
    e.assign_synch_params(&ncgs, &_ncgsnth, &bcgs, &tcgs);
    let nu: Vec<f64> = rshift.iter().map(|g| 3.45e11 / g).collect();
    e.calc_emissivity(&nu);
    // scale as the driver does: j*=g^2/fac, K*=LBH/g
    let fac0: f64 = (0..n).map(|i| e.j[i*4] * (ray.lambda[0] - ray.lambda[i])).sum();
    let lbh = grtrans_core::constants::G * 4e6 * grtrans_core::constants::MSUN / grtrans_core::constants::C2;
    let mut j1 = vec![0.0f64; n];
    let mut k1 = vec![0.0f64; n];
    for i in 0..n {
        j1[i] = e.j[i*4] * rshift[i] * rshift[i] / fac0;
        k1[i] = e.kcoef[i * e.nk] * lbh / rshift[i];
    }
    let jmax = j1.iter().cloned().fold(f64::MIN, f64::max);
    let imax = j1.iter().enumerate().max_by(|a,b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
    let kmax = k1.iter().cloned().fold(f64::MIN, f64::max);
    let ikmax = k1.iter().enumerate().max_by(|a,b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
    println!("  g[0]={:e} fac={fac0:e} LBH={lbh:e}", rshift[0]);
    println!("  j1 max {jmax:e} at {imax}; K1 there {:e}; K1 max {kmax:e} at {ikmax}; j1 there {:e}", k1[imax], j1[ikmax]);
    for i in [imax, ikmax] {
        println!("   idx {i}: j1={:e} K1={:e} j1/K1={:e} jraw={:e} Kraw={:e} bnu_src={:e}", j1[i], k1[i], j1[i]/k1[i], e.j[i*4], e.kcoef[i*e.nk], e.j[i*4]/e.kcoef[i*e.nk]);
    }
    let fac: f64 = (0..n).map(|i| e.j[i*4] * (ray.lambda[0] - ray.lambda[i])).sum();
    let mut sp2 = sp.clone();
    sp2.gmin = vec![1.0; n]; sp2.mu = vec![0.25; n]; sp2.jetalpha = vec![0.02; n];
    let opts = TraceOptions { ename: "POLSYNCHTH".to_string(), iname: "lsoda".to_string(), nvals: 4, freqs: vec![3.45e11], ..Default::default() };
    // replicate the driver integration and print the window profile
    let mut e2 = Emis::select("POLSYNCHTH");
    e2.initialize(n, &rshift, &ang, &cosne);
    e2.assign_synch_params(&ncgs, &_ncgsnth, &bcgs, &tcgs);
    e2.calc_emissivity(&nu);
    let fac2: f64 = (0..n).map(|i| e2.j[i*4] * (ray.lambda[0] - ray.lambda[i])).sum();
    invariant_emis_g2(&mut e2, &rshift);
    for v in e2.j.iter_mut() { *v /= fac2; }
    for v in e2.kcoef.iter_mut() { *v *= lbh; }
    let j1v: Vec<f64> = (0..n).map(|i| e2.j[i*4]).collect();
    let k1v: Vec<f64> = (0..n).map(|i| e2.kcoef[i*e2.nk]).collect();
    let tau: Vec<f64> = grtrans_transfer::calc_opt_depth(&ray.lambda, &k1v).iter().map(|v| -v).collect();
    let (inten, nptsout) = grtrans_transfer::integrate_scalar(grtrans_transfer::Method::Lsoda, &ray.lambda, &j1v, &k1v, &tau, 0.1, 1e-8, 1e-6);
    let nz: Vec<usize> = (0..n).filter(|i| j1v[*i] != 0.0).collect();
    println!("  window: first={} last={} count={} nptsout={}", nz.first().unwrap(), nz.last().unwrap(), nz.len(), nptsout);
    for i in nz.iter().rev().take(4) {
        println!("   idx {i}: lam={:e} j1={:e} I={:e}", ray.lambda[*i], j1v[*i], inten[*i]);
    }
    println!("  I(final)={:e} I_cgs={:e}", inten[nptsout-1], inten[nptsout-1]*fac2*lbh);
    let vals = grtrans::driver::trace_ray(&mut loaded, &ray, status, &sp2, &opts, false);
    println!("  fac={fac:e} I={:?}", &vals[..4]);
}
