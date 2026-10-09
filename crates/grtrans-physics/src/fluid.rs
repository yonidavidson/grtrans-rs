//! Fluid model dispatch and data types.
//!
//! Direct translation of `fluid.f90` (upstream GRTRANS): model constants,
//! the `fluid` / `fluid_args` / `source_params` types, and the per-model
//! `get_*_fluidvars` / `convert_fluidvars_*` dispatchers.
//!
//! Upstream keeps per-model data in module globals (`mbh`, `rin`, ...).
//! The Rust port makes that state explicit in [`LoadedFluid`], which the
//! driver creates once (`load_fluid_model`) and mutates while tracing rays
//! (`get_fluid_vars`).

use crate::models::ffjet::{self, FfjetData};
use crate::models::thindisk::{self, ThindiskState};
use grtrans_core::four_vector::FourVector;

/// Upstream model constants (`fluid.f90` lines 40-43).
pub mod model {
    pub const DUMMY: i32 = 0;
    pub const SPHACC: i32 = 1;
    pub const THINDISK: i32 = 2;
    pub const RIAF: i32 = 3;
    pub const HOTSPOT: i32 = 4;
    pub const PHATDISK: i32 = 5;
    pub const SCHNITTMAN: i32 = 6;
    pub const HARM: i32 = 12;
    pub const FFJET: i32 = 13;
    pub const NUMDISK: i32 = 14;
    pub const THICKDISK: i32 = 15;
    pub const MB09: i32 = 16;
    pub const SARIAF: i32 = 17;
    pub const POWERLAW: i32 = 18;
    pub const HARM3D: i32 = 19;
    pub const HARMPI: i32 = 20;
    pub const TOY: i32 = 21;
    pub const KORAL: i32 = 22;
    pub const KORALNTH: i32 = 23;
    pub const IHARM: i32 = 29;
}

/// Upstream `fluid_args`.
#[derive(Clone, Debug)]
pub struct FluidArgs {
    pub dfile: String,
    pub hfile: String,
    pub gfile: String,
    pub sim: String,
    pub nt: i32,
    pub indf: i32,
    pub nfiles: i32,
    pub jonfix: i32,
    pub nw: i32,
    pub nfreq_tab: i32,
    pub nr: i32,
    pub offset: i32,
    pub dindf: i32,
    pub magcrit: i32,
    pub bl06: i32,
    pub rspot: f64,
    pub r0spot: f64,
    pub n0spot: f64,
    pub tscl: f64,
    pub rscl: f64,
    pub wmin: f64,
    pub wmax: f64,
    pub fmin: f64,
    pub fmax: f64,
    pub rmax: f64,
    pub sigt: f64,
    pub fcol: f64,
    pub mdot: f64,
    pub mbh: f64,
    pub nscl: f64,
    pub nnthscl: f64,
    pub nnthp: f64,
    pub beta: f64,
    pub np: f64,
    pub tp: f64,
    pub rin: f64,
    pub rout: f64,
    pub thin: f64,
    pub thout: f64,
    pub phiin: f64,
    pub phiout: f64,
    pub scalefac: f64,
    pub sigcut: f64,
}

impl Default for FluidArgs {
    fn default() -> Self {
        // defaults matching grtrans_batch.py init()
        FluidArgs {
            dfile: String::new(),
            hfile: String::new(),
            gfile: String::new(),
            sim: String::new(),
            nt: 1,
            indf: 1,
            nfiles: 1,
            jonfix: 1,
            nw: 500,
            nfreq_tab: 100,
            nr: 500,
            offset: 0,
            dindf: 1,
            magcrit: 0,
            bl06: 0,
            rspot: 1.5,
            r0spot: 6.0,
            n0spot: 1e4,
            tscl: 1.0,
            rscl: 6.0,
            wmin: 1e-4,
            wmax: 1e4,
            fmin: 3.33e16,
            fmax: 9e19,
            rmax: 1e4,
            sigt: 0.4,
            fcol: 1.7,
            mdot: 0.1,
            mbh: 10.0,
            nscl: 3e7,
            nnthscl: 8e4,
            nnthp: 2.9,
            beta: 10.0,
            np: 0.0,
            tp: 0.0,
            rin: -1.0,
            rout: 1e8,
            thin: -10.0,
            thout: 10.0,
            phiin: 0.0,
            phiout: 1e4,
            scalefac: 1.0,
            sigcut: 1e10,
        }
    }
}

/// Upstream `source_params`.
#[derive(Clone, Debug, Default)]
pub struct SourceParams {
    pub nfac: f64,
    pub bfac: f64,
    pub mbh: f64,
    pub mdot: f64,
    pub p1: f64,
    pub p2: f64,
    pub gmax: f64,
    pub gminval: f64,
    pub jetalphaval: f64,
    pub muval: f64,
    pub sigcut: f64,
    pub gmin: Vec<f64>,
    pub jetalpha: Vec<f64>,
    pub mu: Vec<f64>,
    pub type_: i32,
}

/// Per-model persistent state (upstream module globals).
#[derive(Clone, Debug, Default)]
pub enum ModelState {
    #[default]
    None,
    Thindisk(ThindiskState),
    Ffjet(FfjetData),
}

/// Result of `load_fluid_model`: the model name and its persistent state.
#[derive(Clone, Debug, Default)]
pub struct LoadedFluid {
    pub name: String,
    pub state: ModelState,
}

/// Upstream `fluid`.
#[derive(Clone, Debug)]
pub struct Fluid {
    pub model: i32,
    pub nfreq: usize,
    pub nrelbin: usize,
    pub rin: f32,
    pub bingammamin: f32,
    pub bingammamax: f32,
    pub sigcut: f32,
    pub gamma: f32,
    /// For THINDISK this stores the temperature T (upstream comment: "rho is
    /// used to store T for thindisk")
    pub rho: Vec<f32>,
    pub p: Vec<f32>,
    pub bmag: Vec<f32>,
    pub rho2: Vec<f32>,
    pub u: Vec<FourVector>,
    pub b: Vec<FourVector>,
    pub npts: usize,
}

impl Fluid {
    pub fn empty() -> Self {
        Fluid {
            model: model::DUMMY,
            nfreq: 0,
            nrelbin: 0,
            rin: 0.0,
            bingammamin: 1.0,
            bingammamax: 1.0,
            sigcut: 1e10,
            gamma: 0.0,
            rho: Vec::new(),
            p: Vec::new(),
            bmag: Vec::new(),
            rho2: Vec::new(),
            u: Vec::new(),
            b: Vec::new(),
            npts: 0,
        }
    }
}

/// Upstream `initialize_fluid_model`.
pub fn initialize_fluid_model(loaded: &LoadedFluid, nup: usize) -> Fluid {
    let mut f = Fluid::empty();
    f.npts = nup;
    f.u = vec![FourVector::flat([0.0; 4]); nup];
    f.b = vec![FourVector::flat([0.0; 4]); nup];
    f.rho = vec![0.0; nup];
    match loaded.name.as_str() {
        "THINDISK" => {
            f.model = model::THINDISK;
        }
        "FFJET" => {
            f.model = model::FFJET;
            f.bmag = vec![0.0; nup];
            f.p = vec![0.0; nup];
        }
        _ => {
            f.model = model::DUMMY;
        }
    }
    f
}

/// Upstream `load_fluid_model` (only the ported models act).
pub fn load_fluid_model(fname: &str, a: f64, args: &FluidArgs) -> LoadedFluid {
    match fname {
        "THINDISK" => {
            let mut state = ThindiskState::default();
            thindisk::init_thindisk(
                &mut state,
                a as f32,
                args.mdot as f32,
                args.mbh as f32,
                args.rin as f32,
                args.rout as f32,
            );
            LoadedFluid {
                name: fname.to_string(),
                state: ModelState::Thindisk(state),
            }
        }
        "FFJET" => {
            let data = ffjet::initialize_ffjet_model(std::path::Path::new(&args.dfile));
            LoadedFluid {
                name: fname.to_string(),
                state: ModelState::Ffjet(data),
            }
        }
        _ => LoadedFluid {
            name: fname.to_string(),
            state: ModelState::None,
        },
    }
}

/// Upstream `get_fluid_vars` (array version): fill fluid variables along a
/// ray.
pub fn get_fluid_vars(
    loaded: &mut LoadedFluid,
    x: &[FourVector],
    k: &[FourVector],
    a: f64,
    f: &mut Fluid,
) {
    match (&mut loaded.state, f.model) {
        (ModelState::Thindisk(state), model::THINDISK) => {
            thindisk::get_thindisk_fluidvars(state, x, k, a as f32, f)
        }
        (ModelState::Ffjet(data), model::FFJET) => {
            let (rho, p, bmag, u, b) = ffjet::ffjet_vals(data, x, a as f32);
            f.rho = rho;
            f.p = p;
            f.bmag = bmag;
            f.u = u;
            f.b = b;
        }
        _ => {
            for v in f.u.iter_mut() {
                *v = FourVector::flat([0.0; 4]);
            }
            for v in f.b.iter_mut() {
                *v = FourVector::flat([0.0; 4]);
            }
        }
    }
}

/// Upstream `convert_fluid_vars` (array version): fluid -> cgs emission
/// quantities `(ncgs, ncgsnth, bcgs, tcgs)`.
pub fn convert_fluid_vars(
    f: &Fluid,
    sp: &SourceParams,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = f.npts;
    match f.model {
        model::THINDISK => {
            // convert_fluidvars_thindisk: tcgs=f%rho, ncgs=1
            let tcgs: Vec<f64> = f.rho.iter().map(|v| *v as f64).collect();
            (vec![1.0; n], vec![0.0; n], vec![0.0; n], tcgs)
        }
        model::FFJET => {
            // convert_fluidvars_ffjet: ncgsnth=rho*nfac, bcgs=bmag*bfac
            let ncgsnth: Vec<f64> = f.rho.iter().map(|v| *v as f64 * sp.nfac).collect();
            let bcgs: Vec<f64> = f.bmag.iter().map(|v| *v as f64 * sp.bfac).collect();
            (vec![0.0; n], ncgsnth, bcgs, vec![0.0; n])
        }
        _ => (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]),
    }
}

/// Upstream `assign_source_params` for `stype='const'`: constant gmin/mu/
/// jetalpha per point.
pub fn assign_source_params(sp: &mut SourceParams) {
    // CONST case (upstream select case(sp%type) case (CONST))
    sp.gmin = vec![sp.gminval; sp.gmin.len().max(1)];
    sp.jetalpha = vec![sp.jetalphaval; sp.jetalpha.len().max(1)];
    sp.mu = vec![sp.muval; sp.mu.len().max(1)];
}
