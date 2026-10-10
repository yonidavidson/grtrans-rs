//! Debug: trace a few SPHACC rays and print intermediate quantities.
use grtrans::driver::{trace_ray, TraceOptions};
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::fluid::{
    convert_fluid_vars, get_fluid_vars, initialize_fluid_model, load_fluid_model, FluidArgs,
    SourceParams,
};

fn main() {
    let args = initialize_pixels(
        true, 1, 0.1, -0.5, 0.0, 0.0025, 1.0, 1.0, 2, 0.0, 400.0, 0.0, 0.0, 10000, 1, 100,
    );
    let fargs = FluidArgs {
        mbh: 1.0,
        ..Default::default()
    };
    let mut loaded = load_fluid_model("SPHACC", 0.0, &fargs);
    let sp = SourceParams {
        mbh: 1.0,
        mdot: 1.0,
        ..Default::default()
    };
    let freqs: Vec<f64> = (0..25)
        .map(|i| 1e8 * ((i as f64) * (1e15f64 / 1e8f64).ln() / 24.0).exp())
        .collect();
    let opts = TraceOptions {
        ename: "SYNCHTHAV".to_string(),
        iname: "lsoda".to_string(),
        nvals: 1,
        freqs,
        ..Default::default()
    };
    if std::env::args().nth(1).as_deref() == Some("dump") {
        let mut ours = vec![0.0f64; 10000 * 25];
        for i in 0..10000 {
            let (ray, status) = initialize_geodesic(&args, i);
            let vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
            for k in 0..25 {
                ours[i * 25 + k] = vals[k];
            }
        }
        let mut bytes = Vec::new();
        for v in &ours {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        std::fs::write("/tmp/sphacc_ours.f64.bin", bytes).unwrap();
        println!("wrote /tmp/sphacc_ours.f64.bin");
        return;
    }
    for pix in [0usize, 5000, 9000] {
        let (ray, status) = initialize_geodesic(&args, pix);
        let n = ray.npts;
        println!(
            "pixel {pix}: status={status} npts={n} r=[{:.4}, {:.4}] u=[{:.6}, {:.6}]",
            ray.x[0].data[1],
            ray.x[n - 1].data[1],
            1.0 / ray.x[0].data[1],
            1.0 / ray.x[n - 1].data[1]
        );
        for k in [0, 1, 2, 5, 20, n / 2, n - 3, n - 1] {
            println!(
                "   k={k:3} r={:.5} u={:.6} lam={:.5e}",
                ray.x[k].data[1],
                1.0 / ray.x[k].data[1],
                ray.lambda[k]
            );
        }
        let mut f = initialize_fluid_model(&loaded, n);
        get_fluid_vars(&mut loaded, &ray.x, &ray.k, 0.0, &mut f);
        for k in [73usize, 74, 75, 76, 77] {
            println!(
                "   raw k={k} r={:.5} u=1/r={:.6} ur={:e} n={:e}",
                ray.x[k].data[1],
                1.0 / ray.x[k].data[1],
                f.u[k].data[1],
                f.rho[k]
            );
        }
        let (ncgs, _ncgsnth, bcgs, tcgs) = convert_fluid_vars(&f, &sp);
        println!(
            "  n[0]={:e} n[last]={:e} B[0]={:e} T[0]={:e} ur[0]={:e}",
            ncgs[0],
            ncgs[n - 1],
            bcgs[0],
            tcgs[0],
            f.u[0].data[1]
        );
        println!(
            "  n mid={:e} B mid={:e} T mid={:e}",
            ncgs[n / 2],
            bcgs[n / 2],
            tcgs[n / 2]
        );
        // instrument the k=0 chain
        {
            use grtrans_core::kerr::comoving_ortho;
            use grtrans_physics::emissivity::{invariant_emis_g2, Emis};
            let mut u = f.u.clone();
            let mut b = f.b.clone();
            let mut kv = ray.k.clone();
            let mut rshift = vec![0.0f64; n];
            let mut ang = vec![0.0f64; n];
            let mut cosne = vec![0.0f64; n];
            for i in 0..n {
                let co = comoving_ortho(
                    ray.x[i].data[1],
                    ray.x[i].data[2],
                    0.0,
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
            }
            let mut e = Emis::select("SYNCHTHAV");
            e.initialize(n, &rshift, &ang, &cosne);
            e.assign_synch_params(&ncgs, &_ncgsnth, &bcgs, &tcgs);
            let nu0: f64 = std::env::args()
                .nth(2)
                .and_then(|v| v.parse().ok())
                .unwrap_or(1e8);
            let nu: Vec<f64> = rshift.iter().map(|g| nu0 / g).collect();
            e.calc_emissivity(&nu);
            let fac: f64 = (0..n)
                .map(|i| e.j[i] * (ray.lambda[0] - ray.lambda[i]))
                .sum();
            println!(
                "  g[0]={:e} g[last]={:e} nu_em[0]={:e}",
                rshift[0],
                rshift[n - 1],
                nu[0]
            );
            println!(
                "  j[0]={:e} j[last]={:e} K[0]={:e} K[last]={:e} fac={:e}",
                e.j[0],
                e.j[n - 1],
                e.kcoef[0],
                e.kcoef[n - 1],
                fac
            );
            invariant_emis_g2(&mut e, &rshift);
            for v in e.j.iter_mut() {
                *v /= fac;
            }
            let lbh = grtrans_core::constants::G * 1.0 * grtrans_core::constants::MSUN
                / grtrans_core::constants::C2;
            for v in e.kcoef.iter_mut() {
                *v *= lbh;
            }
            let alpha_i: Vec<f64> = (0..n).map(|i| e.kcoef[i]).collect();
            let tau: Vec<f64> = grtrans_transfer::calc_opt_depth(&ray.lambda, &alpha_i)
                .iter()
                .map(|v| -v)
                .collect();
            println!(
                "  tau[last]={:e} K1[last]={:e} lbh={:e}",
                tau[n - 1],
                alpha_i[n - 1],
                lbh
            );
            let j1: Vec<f64> = (0..n).map(|i| e.j[i]).collect();
            let inten =
                grtrans_transfer::radtrans_integrate_quadrature(&ray.lambda, &j1, &alpha_i, &tau);
            println!(
                "  nu0={:e} I_quad*fac*lbh = {:e}",
                nu0,
                inten[n - 1] * fac * lbh
            );
            for k in [75usize] {
                println!("   u[{k}]={:?} k[{k}]={:?}", u[k].data, kv[k].data);
            }
            for k in [0usize, 10, 25, 50, 75, 99] {
                println!(
                    "   jscaled[{k}]={:e} K1[{k}]={:e} n={:e} T={:e} B={:e} g={:e} jraw={:e}",
                    e.j[k],
                    alpha_i[k],
                    ncgs[k],
                    tcgs[k],
                    bcgs[k],
                    rshift[k],
                    e.j[k] * fac / (rshift[k] * rshift[k])
                );
            }
        }
        let vals = trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
        println!(
            "  I(k=0,12,24) = {:e} {:e} {:e}",
            vals[0], vals[12], vals[24]
        );
    }
}
