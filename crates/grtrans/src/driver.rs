//! Per-ray driver: fluid, emissivity, transfer, output.
//!
//! Direct translation of `grtrans_driver.f90` (upstream GRTRANS). The
//! integration path (npts > 1) uses `grtrans-transfer`; the single-point
//! path (thin disk, hot spots) computes the intensity directly with the
//! g^3 invariant scaling.

use grtrans_core::constants::{C2, G, MSUN};
use grtrans_core::kerr::comoving_ortho;
use grtrans_geodesics::rays::Ray;
use grtrans_physics::emissivity::{invariant_emis, rotate_emis, Emis, EmisParams};
use grtrans_physics::fluid::{
    convert_fluid_vars, get_fluid_vars, initialize_fluid_model, LoadedFluid,
};

/// Upstream `rad_trans` object (per-ray integration state).
pub struct RadTrans {
    pub iflag: i32,
    pub neq: usize,
    pub npts: usize,
    /// intensity array, row-major (neq x npts)
    pub i: Vec<f64>,
}

impl RadTrans {
    pub fn new(iname: &str, npts: usize, neq: usize) -> Self {
        let iflag = match iname {
            "lsoda" => 0,
            "delo" => 1,
            "formal" => 2,
            "lsodasph" => 3,
            _ => panic!("iname not recognized: {iname}"),
        };
        RadTrans {
            iflag,
            neq,
            npts,
            i: vec![0.0; neq * npts],
        }
    }
}

/// Per-ray result (upstream `ray_set` pixel values).
#[derive(Clone, Debug, Default)]
pub struct RayResult {
    pub values: Vec<f64>,
}

/// Parameters of one driver invocation (upstream `source_params` subset).
#[derive(Clone, Debug)]
pub struct SourceParams {
    pub mdot: f64,
    pub mbh: f64,
}

/// Options for tracing one ray.
#[derive(Clone, Debug)]
pub struct TraceOptions {
    pub ename: String,
    pub iname: String,
    pub nvals: usize,
    pub freqs: Vec<f64>,
}

/// Trace one ray and return the Stokes vector (length `nvals`) for each
/// frequency, flattened as `[freq][stokes]`.
#[allow(clippy::too_many_arguments)]
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
        let out = comoving_ortho(
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
        s2xi[i] = out.s2xi;
        c2xi[i] = out.c2xi;
        ang[i] = out.ang;
        rshift[i] = out.g;
        cosne[i] = out.cosne;
    }
    let mut r = RadTrans::new(&opts.iname, npts, nvals);
    let _ = extra;
    e.initialize(npts, &rshift, &ang, &cosne);
    let (ncgs, bcgs, tcgs) = convert_fluid_vars(&f);
    let _ = (&ncgs, &bcgs);
    e.assign_params(&tcgs);
    let ep = EmisParams::default();
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
                panic!("integrated rays require grtrans-transfer (not yet wired)");
            } else {
                invariant_emis(&mut e, &rshift, 3);
                // grtrans_compute_intensity: r%I = transpose(e%j)
                for i in 0..npts {
                    for q in 0..e.neq {
                        r.i[q * npts + i] = e.j[i * e.neq + q];
                    }
                }
                let _ = &ep;
                let _ = lbh;
            }
            for q in 0..nvals.min(e.neq) {
                out[k * nvals + q] = r.i[q * npts + npts - 1];
            }
        } else {
            // zero intensity
        }
    }
    out
}
