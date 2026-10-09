//! Python bindings for grtrans-rs (PyO3).
//!
//! The public API mirrors the spirit of the upstream `grtrans_batch.py` /
//! `pgrtrans` interface: a single `run_image(**kwargs)` entry point taking
//! the same parameter names as `grtrans_inputs`, plus a low-level
//! `geodesic_ray` helper.
//!
//! Build with maturin:
//! ```bash
//! maturin develop -m crates/grtrans-python/Cargo.toml
//! ```

use numpy::{IntoPyArray, PyArray2, PyArray3, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// grtrans-rs version string.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

fn f64_kw(kw: &Bound<'_, PyDict>, key: &str, default: f64) -> PyResult<f64> {
    match kw.get_item(key)? {
        Some(v) => v.extract::<f64>(),
        None => Ok(default),
    }
}

fn i64_kw(kw: &Bound<'_, PyDict>, key: &str, default: i64) -> PyResult<i64> {
    match kw.get_item(key)? {
        Some(v) => v.extract::<i64>(),
        None => Ok(default),
    }
}

fn str_kw(kw: &Bound<'_, PyDict>, key: &str, default: &str) -> PyResult<String> {
    match kw.get_item(key)? {
        Some(v) => v.extract::<String>(),
        None => Ok(default.to_string()),
    }
}

fn bool_kw(kw: &Bound<'_, PyDict>, key: &str, default: bool) -> PyResult<bool> {
    match kw.get_item(key)? {
        Some(v) => v.extract::<bool>(),
        None => Ok(default),
    }
}

fn vec3_kw(kw: &Bound<'_, PyDict>, key: &str, default: [i64; 3]) -> PyResult<[i64; 3]> {
    match kw.get_item(key)? {
        Some(v) => {
            let v: Vec<i64> = v.extract()?;
            if v.len() != 3 {
                return Err(PyValueError::new_err(format!("{key} must have 3 values")));
            }
            Ok([v[0], v[1], v[2]])
        }
        None => Ok(default),
    }
}

fn vec4_kw(kw: &Bound<'_, PyDict>, key: &str, default: [f64; 4]) -> PyResult<[f64; 4]> {
    match kw.get_item(key)? {
        Some(v) => {
            let v: Vec<f64> = v.extract()?;
            if v.len() != 4 {
                return Err(PyValueError::new_err(format!("{key} must have 4 values")));
            }
            Ok([v[0], v[1], v[2], v[3]])
        }
        None => Ok(default),
    }
}

/// Build the input deck from keyword arguments (same names as
/// `grtrans_batch.grtrans_inputs`).
fn inputs_from_kwargs(kw: &Bound<'_, PyDict>) -> PyResult<grtrans::inputs::Inputs> {
    let mut inp = grtrans::inputs::Inputs::default();
    inp.standard = i64_kw(kw, "standard", inp.standard as i64)? as i32;
    inp.mumin = f64_kw(kw, "mumin", inp.mumin)?;
    inp.mumax = f64_kw(kw, "mumax", inp.mumax)?;
    inp.nmu = i64_kw(kw, "nmu", inp.nmu)?;
    inp.phi0 = f64_kw(kw, "phi0", inp.phi0)?;
    inp.spin = f64_kw(kw, "spin", inp.spin)?;
    inp.uout = f64_kw(kw, "uout", inp.uout)?;
    inp.uin = f64_kw(kw, "uin", inp.uin)?;
    inp.rcut = f64_kw(kw, "rcut", inp.rcut)?;
    inp.nrotype = i64_kw(kw, "nrotype", inp.nrotype as i64)? as i32;
    inp.gridvals = vec4_kw(kw, "gridvals", inp.gridvals)?;
    inp.nn = vec3_kw(kw, "nn", inp.nn)?;
    inp.i1 = i64_kw(kw, "i1", inp.i1)?;
    inp.i2 = i64_kw(kw, "i2", inp.i2)?;
    inp.extra = i64_kw(kw, "extra", inp.extra as i64)? as i32;
    inp.debug = i64_kw(kw, "debug", inp.debug as i64)? as i32;
    inp.fname = str_kw(kw, "fname", &inp.fname)?;
    inp.nmdot = i64_kw(kw, "nmdot", inp.nmdot)?;
    inp.mdotmin = f64_kw(kw, "mdotmin", inp.mdotmin)?;
    inp.mdotmax = f64_kw(kw, "mdotmax", inp.mdotmax)?;
    inp.ename = str_kw(kw, "ename", &inp.ename)?;
    inp.mbh = f64_kw(kw, "mbh", inp.mbh)?;
    inp.nfreq = i64_kw(kw, "nfreq", inp.nfreq)?;
    inp.fmin = f64_kw(kw, "fmin", inp.fmin)?;
    inp.fmax = f64_kw(kw, "fmax", inp.fmax)?;
    inp.muval = f64_kw(kw, "muval", inp.muval)?;
    inp.gmin = f64_kw(kw, "gmin", inp.gmin)?;
    inp.gmax = f64_kw(kw, "gmax", inp.gmax)?;
    inp.p1 = f64_kw(kw, "p1", inp.p1)?;
    inp.p2 = f64_kw(kw, "p2", inp.p2)?;
    inp.jetalpha = f64_kw(kw, "jetalpha", inp.jetalpha)?;
    inp.stype = str_kw(kw, "stype", &inp.stype)?;
    inp.use_geokerr = bool_kw(kw, "use_geokerr", inp.use_geokerr)?;
    inp.nvals = i64_kw(kw, "nvals", inp.nvals)?;
    inp.iname = str_kw(kw, "iname", &inp.iname)?;
    inp.cflag = i64_kw(kw, "cflag", inp.cflag as i64)? as i32;
    inp.fdfile = str_kw(kw, "jdfile", &inp.fdfile)?;
    inp.fdfile = str_kw(kw, "dfile", &inp.fdfile)?;
    inp.fhfile = str_kw(kw, "hfile", &inp.fhfile)?;
    inp.fhfile = str_kw(kw, "hhfile", &inp.fhfile)?;
    inp.fsim = str_kw(kw, "fsim", &inp.fsim)?;
    inp.fmdot = f64_kw(kw, "tmdot", inp.fmdot)?;
    inp.ftscl = f64_kw(kw, "ntscl", inp.ftscl)?;
    inp.frscl = f64_kw(kw, "nrscl", inp.frscl)?;
    inp.frin = f64_kw(kw, "frin", inp.frin)?;
    inp.frout = f64_kw(kw, "frout", inp.frout)?;
    inp.sigcut = f64_kw(kw, "sigcut", inp.sigcut)?;
    Ok(inp)
}

/// Run a GRTRANS image calculation.
///
/// Keyword arguments match `grtrans_batch.grtrans_inputs` (e.g. `fname`,
/// `standard`, `spin`, `nn`, `gridvals`, `ename`, `nvals`, `nfreq`, `fmin`,
/// `fmax`, `mumin`, `mumax`, `nmu`, `uout`, `mbh`, `tmdot`, `jdfile`).
///
/// Returns `(ab, ivals, nu)`:
/// * `ab`: float64 array (npix, 2) pixel locations
/// * `ivals`: float64 array (npix, nvals, nimages)
/// * `nu`: float64 array (nimages,) observed frequencies
#[pyfunction]
#[pyo3(signature = (**kwargs))]
fn run_image<'py>(
    py: Python<'py>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<(
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray3<f64>>,
    Bound<'py, numpy::PyArray1<f64>>,
)> {
    let empty = PyDict::new(py);
    let kw = kwargs.unwrap_or(&empty);
    let inp = inputs_from_kwargs(kw)?;
    let images = py.detach(|| grtrans::run::run(&inp));
    if images.is_empty() {
        return Err(PyValueError::new_err("no images produced"));
    }
    let npix = images[0].nx * images[0].ny;
    let nvals = images[0].nvals;
    let nimg = images.len();
    let mut ab = vec![0.0f64; npix * 2];
    let mut ivals = vec![0.0f64; npix * nvals * nimg];
    let mut nu = vec![0.0f64; nimg];
    for (k, img) in images.iter().enumerate() {
        nu[k] = img.nu;
        for i in 0..npix {
            ab[i * 2] = img.ab[i][0] as f64;
            ab[i * 2 + 1] = img.ab[i][1] as f64;
            for q in 0..nvals {
                ivals[(i * nvals + q) * nimg + k] = img.ivals[i * nvals + q] as f64;
            }
        }
    }
    let ab = ab.into_pyarray(py).reshape([npix, 2])?;
    let ivals = ivals.into_pyarray(py).reshape([npix, nvals, nimg])?;
    let nu = nu.into_pyarray(py);
    Ok((ab, ivals, nu))
}

/// Trace a single geodesic ray (low-level helper for tests and debugging).
///
/// Returns `(x, k, lambda)` with `x` the BL coordinates (npts, 4),
/// `k` the covariant wave vector (npts, 4), `lambda` the affine parameter.
#[pyfunction]
#[pyo3(signature = (alpha, beta, mu0, spin, standard=1, uout=0.01, nup=400))]
fn geodesic_ray<'py>(
    py: Python<'py>,
    alpha: f64,
    beta: f64,
    mu0: f64,
    spin: f64,
    standard: i32,
    uout: f64,
    nup: usize,
) -> PyResult<(
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, numpy::PyArray1<f64>>,
)> {
    let args = grtrans_geodesics::rays::initialize_pixels(
        true, standard, mu0, -0.5, spin, uout, 1.0, 1.0, 2, alpha, alpha, beta, beta, 1, 1, nup,
    );
    let (ray, _status) = py.detach(|| grtrans_geodesics::rays::initialize_geodesic(&args, 0));
    let n = ray.npts;
    let mut x = vec![0.0f64; n * 4];
    let mut k = vec![0.0f64; n * 4];
    let mut lambda = vec![0.0f64; n];
    for i in 0..n {
        for c in 0..4 {
            x[i * 4 + c] = ray.x[i].data[c];
            k[i * 4 + c] = ray.k[i].data[c];
        }
        lambda[i] = ray.lambda[i];
    }
    let x = x.into_pyarray(py).reshape([n, 4])?;
    let k = k.into_pyarray(py).reshape([n, 4])?;
    let lambda = lambda.into_pyarray(py);
    Ok((x, k, lambda))
}

#[pymodule]
fn grtrans_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(run_image, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_ray, m)?)?;
    Ok(())
}
