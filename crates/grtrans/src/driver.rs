//! Per-ray driver: fluid, emissivity, transfer, output.
//!
//! Direct translation of `grtrans_driver.f90` (upstream GRTRANS).

use grtrans_core::constants::{C2, G, MSUN};
use grtrans_core::kerr::comoving_ortho;
use grtrans_geodesics::rays::Ray;
use grtrans_physics::emissivity::{
    invariant_emis, invariant_emis_g2, rotate_emis, Emis, EmisParams,
};
pub use grtrans_physics::fluid::SourceParams;
use grtrans_physics::fluid::{
    convert_fluid_vars, get_fluid_vars, initialize_fluid_model, LoadedFluid,
};
use grtrans_transfer::{calc_opt_depth, integrate, Method};

/// Upstream `rad_trans` object (per-ray integration state).
pub struct RadTrans {
    pub method: Method,
    pub neq: usize,
    pub npts: usize,
    /// intensity array, row-major (neq x npts)
    pub i: Vec<f64>,
}

impl RadTrans {
    pub fn new(iname: &str, npts: usize, neq: usize) -> Self {
        RadTrans {
            method: Method::from_name(iname),
            neq,
            npts,
            i: vec![0.0; neq * npts],
        }
    }
}

/// Options for tracing one ray.
#[derive(Clone, Debug)]
pub struct TraceOptions {
    pub ename: String,
    pub iname: String,
    pub nvals: usize,
    pub freqs: Vec<f64>,
    /// upstream `thin` (delo threshold)
    pub thin: f64,
    /// upstream integration tolerances
    pub hmax: f64,
    pub oatol: f64,
    pub ortol: f64,
}

impl Default for TraceOptions {
    fn default() -> Self {
        TraceOptions {
            ename: "BB".to_string(),
            iname: "lsoda".to_string(),
            nvals: 1,
            freqs: Vec::new(),
            thin: 1e-2,
            hmax: 0.1,
            oatol: 1e-8,
            ortol: 1e-6,
        }
    }
}

/// Trace one ray and return the Stokes vector (length `nvals`) per
/// frequency, flattened as `[freq][stokes]`.
pub fn trace_ray(
    loaded: &mut LoadedFluid,
    ray: &Ray,
    status: i32,
    sp: &SourceParams,
    opts: &TraceOptions,
    extra: bool,
) -> Vec<f64> {
    let nfreq = opts.freqs.len();
    let nvals = opts.nvals;
    let mut out = vec![0.0f64; nfreq * nvals];
    if status != 1 || ray.npts == 0 {
        return out;
    }
    let npts = ray.npts;
    let a = ray.a;
    let mut f = initialize_fluid_model(loaded, npts);
    let mut e = Emis::select(&opts.ename);
    get_fluid_vars(loaded, &ray.x, &ray.k, a, &mut f);
    // comoving orthonormal frame
    let mut u = f.u.clone();
    let mut b = f.b.clone();
    let mut kv = ray.k.clone();
    let mut s2xi = vec![0.0; npts];
    let mut c2xi = vec![0.0; npts];
    let mut ang = vec![0.0; npts];
    let mut rshift = vec![0.0; npts];
    let mut cosne = vec![0.0; npts];
    for i in 0..npts {
        let co = comoving_ortho(
            ray.x[i].data[1],
            ray.x[i].data[2],
            a,
            ray.alpha,
            ray.beta,
            ray.mu0,
            &mut u[i],
            &mut b[i],
            &mut kv[i],
        );
        s2xi[i] = co.s2xi;
        c2xi[i] = co.c2xi;
        ang[i] = co.ang;
        rshift[i] = co.g;
        cosne[i] = co.cosne;
    }
    let mut r = RadTrans::new(&opts.iname, npts, nvals);
    let _ = extra;
    e.initialize(npts, &rshift, &ang, &cosne);
    let (ncgs, ncgsnth, bcgs, tcgs) = convert_fluid_vars(&f, sp);
    // assign_emis_params (all ported types)
    e.assign_params(&tcgs);
    e.assign_synch_params(&ncgs, &ncgsnth, &bcgs, &tcgs);
    // emis_model for the power-law synchrotron
    if e.type_ == grtrans_physics::emissivity::etype::EPOLSYNCHPL
        || e.type_ == grtrans_physics::emissivity::etype::ESYNCHPL
    {
        let mut gmin = sp.gmin.clone();
        if gmin.len() < npts {
            gmin = vec![sp.gminval; npts];
        }
        e.set_synchpl_model(sp.p1, &gmin, sp.gmax);
    }
    let ep = EmisParams::default();
    let _ = &ep;
    let lbh = G * sp.mbh * MSUN / C2;
    for k in 0..nfreq {
        let nu: Vec<f64> = rshift.iter().map(|g| opts.freqs[k] / g).collect();
        e.calc_emissivity(&nu);
        let fac = if npts > 1 {
            let mut acc = 0.0;
            for i in 0..npts {
                acc += e.j[i * e.neq] * (ray.lambda[0] - ray.lambda[i]);
            }
            acc
        } else {
            e.j[0] / npts as f64
        };
        if fac / 1e16 > 0.0 {
            // upstream rotates the emission/absorption frame for BOTH the
            // single-point and integrated paths (grtrans_driver.f90 line 201)
            if e.neq == 4 {
                rotate_emis(&mut e, &s2xi, &c2xi);
            }
            if npts != 1 {
                invariant_emis_g2(&mut e, &rshift);
                // scale emission close to 1 and put anu in cgs units
                for v in e.j.iter_mut() {
                    *v /= fac;
                }
                for v in e.kcoef.iter_mut() {
                    *v *= lbh;
                }
                let alpha_i: Vec<f64> = (0..npts).map(|i| e.kcoef[i * e.nk]).collect();
                let tau0 = calc_opt_depth(&ray.lambda, &alpha_i);
                let tau: Vec<f64> = tau0.iter().map(|v| -v).collect();
                if e.neq == 4 {
                    let jv: Vec<[f64; 4]> = (0..npts)
                        .map(|i| [e.j[i * 4], e.j[i * 4 + 1], e.j[i * 4 + 2], e.j[i * 4 + 3]])
                        .collect();
                    let kv: Vec<[f64; 7]> = (0..npts)
                        .map(|i| {
                            let b = i * e.nk;
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
                    let (intensity, nptsout) = integrate(
                        r.method,
                        &ray.lambda,
                        &jv,
                        &kv,
                        &tau,
                        opts.thin,
                        opts.hmax,
                        opts.oatol,
                        opts.ortol,
                    );
                    r.npts = nptsout;
                    for i in 0..npts {
                        for q in 0..4 {
                            r.i[q * npts + i] = intensity[i][q] * fac * lbh;
                        }
                    }
                } else {
                    let j1: Vec<f64> = (0..npts).map(|i| e.j[i]).collect();
                    let k1: Vec<f64> = (0..npts).map(|i| e.kcoef[i * e.nk]).collect();
                    // upstream: lsoda solves the scalar ODE over the trimmed
                    // window; delo/formal use the trapezoidal quadrature
                    let (intensity, nptsout) = grtrans_transfer::integrate_scalar(
                        r.method,
                        &ray.lambda,
                        &j1,
                        &k1,
                        &tau,
                        opts.hmax,
                        opts.oatol,
                        opts.ortol,
                    );
                    r.npts = nptsout;
                    for i in 0..npts {
                        r.i[i] = intensity[i] * fac * lbh;
                    }
                }
            } else {
                invariant_emis(&mut e, &rshift, 3);
                // grtrans_compute_intensity: r%I = transpose(e%j)
                for i in 0..npts {
                    for q in 0..e.neq {
                        r.i[q * npts + i] = e.j[i * e.neq + q];
                    }
                }
            }
            for q in 0..nvals.min(e.neq) {
                out[k * nvals + q] = r.i[q * npts + r.npts - 1];
            }
        }
    }
    out
}
