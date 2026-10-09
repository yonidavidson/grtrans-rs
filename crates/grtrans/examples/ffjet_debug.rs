//! Debug: trace one FFJET ray and print intermediate quantities.
use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_core::kerr::comoving_ortho;
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{
    convert_fluid_vars, get_fluid_vars, initialize_fluid_model, load_fluid_model, FluidArgs,
    SourceParams,
};
use grtrans_physics::polsynch::polsynchpl;

fn main() {
    let spin = 0.998;
    let mu0 = 0.906;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -40.0, 20.0, -20.0, 40.0, 100, 100, 400,
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
    let sp = SourceParams {
        nfac: 2.0,
        bfac: 70.0,
        mbh: 3.4e9,
        mdot: 1.5e15,
        p1: 3.5,
        gmax: 1e5,
        gminval: 100.0,
        jetalphaval: 0.02,
        muval: 0.25,
        ..Default::default()
    };
    for pix in [0usize, 2500, 5000, 7500] {
        let (ray, status) = initialize_geodesic(&args, pix);
        println!(
            "pixel {pix}: status={status} npts={} r0={} r_last={}",
            ray.npts,
            ray.x[0].data[1],
            ray.x[ray.npts - 1].data[1]
        );
        let mut f = initialize_fluid_model(&loaded, ray.npts);
        get_fluid_vars(&mut loaded, &ray.x, &ray.k, spin, &mut f);
        let (ncgs, ncgsnth, bcgs, tcgs) = convert_fluid_vars(&f, &sp);
        let rho_max = f.rho.iter().cloned().fold(f32::MIN, f32::max);
        let b_max = f.bmag.iter().cloned().fold(f32::MIN, f32::max);
        println!(
            "  rho_max={rho_max:e} bmag_max={b_max:e} bcgs_max={:e} ncgsnth_max={:e}",
            bcgs.iter().cloned().fold(f64::MIN, f64::max),
            ncgsnth.iter().cloned().fold(f64::MIN, f64::max)
        );
        let _ = (&ncgs, &tcgs);
        // comoving frame + emissivity
        {
            let mut u = f.u.clone();
            let mut b = f.b.clone();
            let mut kv = ray.k.clone();
            let n = ray.npts;
            let mut rshift = vec![0.0f64; n];
            let mut ang = vec![0.0f64; n];
            let mut cosne = vec![0.0f64; n];
            let mut s2xi_dbg = vec![0.0f64; n];
            let mut c2xi_dbg = vec![0.0f64; n];
            for i in 0..n {
                let co = comoving_ortho(
                    ray.x[i].data[1],
                    ray.x[i].data[2],
                    spin,
                    ray.alpha,
                    ray.beta,
                    ray.mu0,
                    &mut u[i],
                    &mut b[i],
                    &mut kv[i],
                );
                rshift[i] = co.g;
                ang[i] = co.ang;
                cosne[i] = co.cosne;
                s2xi_dbg[i] = co.s2xi;
                c2xi_dbg[i] = co.c2xi;
            }
            let nu: Vec<f64> = rshift.iter().map(|g| 3.45e11 / g).collect();
            let gmin = vec![100.0f64; n];
            let e = polsynchpl(&nu, &ncgsnth, &bcgs, &ang, &vec![3.5f64; n], &gmin, 1e5);
            println!(
                "  g range: {:e}..{:e}",
                rshift.iter().cloned().fold(f64::MAX, f64::min),
                rshift.iter().cloned().fold(f64::MIN, f64::max)
            );
            println!(
                "  nu[0]={:e} ang[0]={:e} cosne[0]={:e}",
                nu[0], ang[0], cosne[0]
            );
            println!("  emis[0]={:?}", e[0]);
            let jmax = e.iter().map(|r| r[0]).fold(f64::MIN, f64::max);
            println!("  jI max={jmax:e}");
            let _ = cosne;
            // replicate driver integration steps
            let n = ray.npts;
            let mut e2 = grtrans_physics::emissivity::Emis::select("POLSYNCHPL");
            e2.initialize(n, &rshift, &ang, &cosne);
            e2.assign_synch_params(&ncgs, &ncgsnth, &bcgs, &tcgs);
            e2.set_synchpl_model(3.5, &vec![100.0f64; n], 1e5);
            e2.calc_emissivity(&nu);
            let fac: f64 = (0..n)
                .map(|i| e2.j[i * 4] * (ray.lambda[0] - ray.lambda[i]))
                .sum();
            println!("  fac={fac:e}");
            if e2.neq == 4 {
                grtrans_physics::emissivity::rotate_emis(&mut e2, &s2xi_dbg, &c2xi_dbg);
            }
            grtrans_physics::emissivity::invariant_emis_g2(&mut e2, &rshift);
            for v in e2.j.iter_mut() {
                *v /= fac;
            }
            let lbh = grtrans_core::constants::G * 3.4e9 * grtrans_core::constants::MSUN
                / grtrans_core::constants::C2;
            for v in e2.kcoef.iter_mut() {
                *v *= lbh;
            }
            let alpha_i: Vec<f64> = (0..n).map(|i| e2.kcoef[i * e2.nk]).collect();
            let tau0 = grtrans_transfer::calc_opt_depth(&ray.lambda, &alpha_i);
            let tau: Vec<f64> = tau0.iter().map(|v| -v).collect();
            let jv: Vec<[f64; 4]> = (0..n)
                .map(|i| {
                    [
                        e2.j[i * 4],
                        e2.j[i * 4 + 1],
                        e2.j[i * 4 + 2],
                        e2.j[i * 4 + 3],
                    ]
                })
                .collect();
            let kv: Vec<[f64; 7]> = (0..n)
                .map(|i| {
                    let b = i * e2.nk;
                    [
                        e2.kcoef[b],
                        e2.kcoef[b + 1],
                        e2.kcoef[b + 2],
                        e2.kcoef[b + 3],
                        e2.kcoef[b + 4],
                        e2.kcoef[b + 5],
                        e2.kcoef[b + 6],
                    ]
                })
                .collect();
            let (intensity, nptsout) = grtrans_transfer::integrate(
                grtrans_transfer::Method::Lsoda,
                &ray.lambda,
                &jv,
                &kv,
                &tau,
                1e-2,
                0.1,
                1e-8,
                1e-6,
            );
            println!("  nptsout={nptsout} I(last)={:?}", intensity[nptsout - 1]);
            println!(
                "  lambda[0]={:e} lambda[last]={:e} tau[last]={:e} K0[last]={:e} j0[last]={:e}",
                ray.lambda[0],
                ray.lambda[n - 1],
                tau[n - 1],
                kv[n - 1][0],
                jv[n - 1][0]
            );
        }
        // full trace
        let opts = TraceOptions {
            ename: "POLSYNCHPL".to_string(),
            iname: "lsoda".to_string(),
            nvals: 4,
            freqs: vec![3.45e11],
            ..Default::default()
        };
        let mut sp2 = sp.clone();
        sp2.gmin = vec![100.0; ray.npts];
        sp2.mu = vec![0.25; ray.npts];
        sp2.jetalpha = vec![0.02; ray.npts];
        let vals = trace_ray(&mut loaded, &ray, status, &sp2, &opts, false);
        println!("  trace: {vals:?}");
    }
}
