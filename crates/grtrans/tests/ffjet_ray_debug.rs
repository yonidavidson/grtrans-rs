//! Per-ray FFJET debug comparison (upstream `debug=1` geodebug output).
//!
//! Compares the full chain for pixel 6035 (1-based): fluid variables,
//! comoving frame, emissivity (after rotation/invariant/scaling), optical
//! depth and final intensity.

use grtrans::driver::TraceOptions;
use grtrans_core::kerr::comoving_ortho;
use grtrans_geodesics::rays::{initialize_geodesic, initialize_pixels};
use grtrans_physics::emissivity::{invariant_emis_g2, rotate_emis, Emis};
use grtrans_physics::fluid::{
    convert_fluid_vars, get_fluid_vars, initialize_fluid_model, load_fluid_model, FluidArgs,
    SourceParams,
};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/ffjet_pixel6035_geodebug.txt"
);

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn ffjet_ray_debug_chain_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("geodebug fixture missing");
    let lines: Vec<Vec<f64>> = text
        .lines()
        .map(|l| l.split_whitespace().map(|t| t.parse().unwrap()).collect())
        .collect();
    let hdr = &lines[0];
    let npts = hdr[0] as usize;
    let fac_f = hdr[8];
    let mu0 = hdr[10];
    assert_eq!(npts, 454);

    // reconstruct the pixel from the header alpha/beta
    let alpha_f = hdr[2];
    let beta_f = hdr[3];

    let spin = 0.998;
    let args = initialize_pixels(
        true, 1, mu0, -0.5, spin, 0.01, 1.0, 1.0, 2, -40.0, 20.0, -20.0, 40.0, 100, 100, 400,
    );
    // find the pixel with matching alpha/beta
    let mut pix = None;
    for (i, p) in args.pixels.iter().enumerate() {
        if (p.alpha - alpha_f).abs() < 1e-12 && (p.beta - beta_f).abs() < 1e-12 {
            pix = Some(i);
            break;
        }
    }
    let pix = pix.expect("pixel with matching alpha/beta not found");

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
    let (ray, status) = initialize_geodesic(&args, pix);
    assert_eq!(status, 1);

    let mut f = initialize_fluid_model(&loaded, ray.npts);
    get_fluid_vars(&mut loaded, &ray.x, &ray.k, spin, &mut f);
    let (ncgs, ncgsnth, bcgs, tcgs) = convert_fluid_vars(
        &f,
        &SourceParams {
            nfac: 2.0,
            bfac: 70.0,
            ..Default::default()
        },
    );

    // stage 1: ray geometry
    let mut worst = 0.0f64;
    let mut worst_desc = String::new();
    for i in 0..npts {
        let checks: [(&str, f64, f64); 9] = [
            ("lambda", ray.lambda[i], lines[1][i]),
            ("t", ray.x[i].data[0], lines[2][i]),
            ("r", ray.x[i].data[1], lines[3][i]),
            ("theta", ray.x[i].data[2], lines[4][i]),
            ("phi", ray.x[i].data[3], lines[5][i]),
            ("kt", ray.k[i].data[0], lines[6][i]),
            ("kr", ray.k[i].data[1], lines[7][i]),
            ("kth", ray.k[i].data[2], lines[8][i]),
            ("kph", ray.k[i].data[3], lines[9][i]),
        ];
        for (name, a, b) in checks {
            let e = (a - b).abs() / b.abs().max(1e-9);
            if e > worst {
                worst = e;
                worst_desc = format!("{name} i={i}: rust {a:e} fort {b:e}");
            }
        }
    }
    println!("stage ray: worst {worst:e} ({worst_desc})");
    assert!(worst < 1e-7, "ray stage: {worst:e} ({worst_desc})");

    // stage 2: fluid + comoving frame
    let mut u = f.u.clone();
    let mut b = f.b.clone();
    let mut kv = ray.k.clone();
    let n = ray.npts;
    let mut rshift = vec![0.0f64; n];
    let mut ang = vec![0.0f64; n];
    let mut cosne = vec![0.0f64; n];
    let mut s2xi = vec![0.0f64; n];
    let mut c2xi = vec![0.0f64; n];
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
        s2xi[i] = co.s2xi;
        c2xi[i] = co.c2xi;
    }
    let mut worst = 0.0f64;
    let mut worst_desc = String::new();
    let mut name_errors: std::collections::HashMap<String, (f64, usize, f64, f64)> =
        std::collections::HashMap::new();
    for i in 0..n {
        let checks: [(&str, f64, f64); 19] = [
            ("rho", f.rho[i] as f64, lines[12][i]),
            ("p", f.p[i] as f64, lines[13][i]),
            ("bmag", f.bmag[i] as f64, lines[14][i]),
            ("ut", u[i].data[0], lines[15][i]),
            ("ur", u[i].data[1], lines[16][i]),
            ("uth", u[i].data[2], lines[17][i]),
            ("uph", u[i].data[3], lines[18][i]),
            ("bt", b[i].data[0], lines[19][i]),
            ("br", b[i].data[1], lines[20][i]),
            ("bth", b[i].data[2], lines[21][i]),
            ("bph", b[i].data[3], lines[22][i]),
            ("ncgs", ncgs[i], lines[23][i]),
            ("tcgs", tcgs[i], lines[24][i]),
            ("bcgs", bcgs[i], lines[25][i]),
            ("g", rshift[i], lines[26][i]),
            ("incang", ang[i], lines[27][i]),
            ("s2xi", s2xi[i], lines[44][i]),
            ("c2xi", c2xi[i], lines[45][i]),
            (
                "kb_ang",
                grtrans_core::kerr::calc_kb_ang(
                    &kv[i],
                    &b[i],
                    &u[i],
                    ray.x[i].data[1],
                    ray.x[i].data[2],
                    spin,
                ),
                lines[46][i],
            ),
        ];
        for (name, a, bb) in checks {
            // s2xi/c2xi are ill-conditioned when the polarization angle is
            // near zero (division by aadotbp^2+bpdotbb^2); use an absolute
            // tolerance for them.
            let scale = if name == "s2xi" || name == "c2xi" {
                1.0
            } else {
                bb.abs().max(1e-3)
            };
            let e = (a - bb).abs() / scale;
            let entry = name_errors
                .entry(name.to_string())
                .or_insert((0.0, 0usize, 0.0, 0.0));
            if e > entry.0 {
                *entry = (e, i, a, bb);
            }
            if e > worst {
                worst = e;
                worst_desc = format!("{name} i={i}: rust {a:e} fort {bb:e}");
            }
        }
    }
    for (name, (e, i, a, bb)) in &name_errors {
        println!("  {name:8}: worst {e:e} at i={i} (rust {a:e} fort {bb:e})");
    }
    // focused debug at the worst incang point
    {
        let i = 329usize;
        let kht = lines[55][i];
        let khr = lines[56][i];
        let khth = lines[57][i];
        let khph = lines[58][i];
        let bht = lines[59][i];
        let bhr = lines[60][i];
        let bhth = lines[61][i];
        let bhph = lines[62][i];
        println!("i=329 fort khat: {kht:e} {khr:e} {khth:e} {khph:e}");
        println!("i=329 fort bhat: {bht:e} {bhr:e} {bhth:e} {bhph:e}");
        println!(
            "i=329 rust b: {:?} u: {:?} k: {:?}",
            b[i].data, u[i].data, kv[i].data
        );
        // rust frame internals via a fresh call on copies
        let mut u2 = u[i];
        let mut b2 = b[i];
        let mut k2 = kv[i];
        let co = comoving_ortho(
            ray.x[i].data[1],
            ray.x[i].data[2],
            spin,
            ray.alpha,
            ray.beta,
            ray.mu0,
            &mut u2,
            &mut b2,
            &mut k2,
        );
        println!(
            "i=329 rust ang={:e} s2xi={:e} c2xi={:e} g={:e}",
            co.ang, co.s2xi, co.c2xi, co.g
        );
    }
    println!("stage fluid/frame: worst {worst:e} ({worst_desc})");
    assert!(worst < 1e-4, "fluid/frame stage: {worst:e} ({worst_desc})");

    // stage 3: emissivity after rotation + invariant + scaling
    let mut sp = SourceParams {
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
    sp.gmin = vec![100.0; n];
    sp.mu = vec![0.25; n];
    sp.jetalpha = vec![0.02; n];
    let mut e = Emis::select("POLSYNCHPL");
    e.initialize(n, &rshift, &ang, &cosne);
    e.assign_synch_params(&ncgs, &ncgsnth, &bcgs, &tcgs);
    e.set_synchpl_model(3.5, &sp.gmin, 1e5);
    let nu: Vec<f64> = rshift.iter().map(|g| 3.45e11 / g).collect();
    e.calc_emissivity(&nu);
    let fac: f64 = (0..n)
        .map(|i| e.j[i * 4] * (ray.lambda[0] - ray.lambda[i]))
        .sum();
    // NOTE: the upstream geodebug header writes `fac` before it is computed
    // (uninitialized), so it cannot be compared; `fac_f` is unused.
    println!("fac: rust {fac:e} (header value is uninitialized upstream: {fac_f:e})");
    let _ = fac_f;
    rotate_emis(&mut e, &s2xi, &c2xi);
    invariant_emis_g2(&mut e, &rshift);
    for v in e.j.iter_mut() {
        *v /= fac;
    }
    let lbh = grtrans_core::constants::G * 3.4e9 * grtrans_core::constants::MSUN
        / grtrans_core::constants::C2;
    for v in e.kcoef.iter_mut() {
        *v *= lbh;
    }
    let mut worst = 0.0f64;
    let mut worst_desc = String::new();
    for i in 0..n {
        let checks: [(&str, f64, f64); 8] = [
            ("jI", e.j[i * 4], lines[28][i]),
            ("KI", e.kcoef[i * 7], lines[29][i]),
            ("jQ", e.j[i * 4 + 1], lines[30][i]),
            ("jV", e.j[i * 4 + 3], lines[31][i]),
            ("KQ", e.kcoef[i * 7 + 1], lines[32][i]),
            ("KV", e.kcoef[i * 7 + 3], lines[33][i]),
            ("rhoQ", e.kcoef[i * 7 + 4], lines[34][i]),
            ("rhoV", e.kcoef[i * 7 + 6], lines[35][i]),
        ];
        for (name, a, bb) in checks {
            let scale = bb.abs().max(1e-12);
            let err = (a - bb).abs() / scale;
            if err > worst {
                worst = err;
                worst_desc = format!("{name} i={i}: rust {a:e} fort {bb:e}");
            }
        }
    }
    println!("stage emissivity: worst {worst:e} ({worst_desc})");
    assert!(worst < 1e-3, "emissivity stage: {worst:e} ({worst_desc})");

    // stage 4: optical depth + final intensity
    let alpha_i: Vec<f64> = (0..n).map(|i| e.kcoef[i * 7]).collect();
    let tau: Vec<f64> = grtrans_transfer::calc_opt_depth(&ray.lambda, &alpha_i)
        .iter()
        .map(|v| -v)
        .collect();
    let mut worst = 0.0f64;
    let mut worst_desc = String::new();
    for i in 0..n {
        let err = rel(tau[i], lines[40][i]);
        if err > worst {
            worst = err;
            worst_desc = format!("tau i={i}: rust {:e} fort {:e}", tau[i], lines[40][i]);
        }
    }
    println!("stage tau: worst {worst:e} ({worst_desc})");

    // final intensity through the same driver call, scaled back
    let opts = TraceOptions {
        ename: "POLSYNCHPL".to_string(),
        iname: "lsoda".to_string(),
        nvals: 4,
        freqs: vec![3.45e11],
        ..Default::default()
    };
    // compare both of my integrators on the same coefficients
    {
        let jv: Vec<[f64; 4]> = (0..n)
            .map(|i| [e.j[i * 4], e.j[i * 4 + 1], e.j[i * 4 + 2], e.j[i * 4 + 3]])
            .collect();
        let kvv: Vec<[f64; 7]> = (0..n)
            .map(|i| {
                let b = i * 7;
                [
                    e.kcoef[b],
                    e.kcoef[b + 1],
                    e.kcoef[b + 2],
                    e.kcoef[b + 3],
                    e.kcoef[b + 4],
                    e.kcoef[b + 5],
                    e.kcoef[b + 6],
                ]
            })
            .collect();
        let (int_delo, _) = grtrans_transfer::integrate(
            grtrans_transfer::Method::Delo,
            &ray.lambda,
            &jv,
            &kvv,
            &tau,
            1e-2,
            0.1,
            1e-8,
            1e-6,
        );
        let (int_lsoda, _) = grtrans_transfer::integrate(
            grtrans_transfer::Method::Lsoda,
            &ray.lambda,
            &jv,
            &kvv,
            &tau,
            1e-2,
            0.1,
            1e-8,
            1e-6,
        );
        println!(
            "raw delo I0={:e} lsoda I0={:e}",
            int_delo[n - 1][0],
            int_lsoda[n - 1][0]
        );
        // thin-limit trapezoid of j over lambda
        let mut acc = 0.0;
        for i in 1..n {
            acc += 0.5 * (jv[i][0] + jv[i - 1][0]) * (ray.lambda[i] - ray.lambda[i - 1]);
        }
        println!("raw thin-limit integral I0={:e}", acc);
    }
    let vals = grtrans::driver::trace_ray(&mut loaded, &ray, status, &sp, &opts, false);
    let norm = fac * lbh;
    let i_fort = [
        lines[36][n - 1],
        lines[37][n - 1],
        lines[38][n - 1],
        lines[39][n - 1],
    ];
    for q in 0..4 {
        let ours_n = vals[q] / norm;
        let r = rel(ours_n, i_fort[q]);
        println!("I{q}: rust {ours_n:e} fort {:e} rel {r:e}", i_fort[q]);
    }
}
