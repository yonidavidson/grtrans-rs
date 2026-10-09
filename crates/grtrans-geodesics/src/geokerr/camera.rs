//! INITIALIZE_CAMERA_GEOKERR: camera pixels to impact parameters and
//! constants of motion. Translation of `geokerr_wrapper.f` lines 1-354.
//!
//! GRTRANS only uses `standard = 1` (radial tracing) and `standard = 2`
//! (polar-angle tracing to the equatorial plane); the manual/interactive
//! `standard = 0/99` paths of the original geokerr program are not
//! reachable from GRTRANS and are not ported (documented in
//! docs/PORTING_MATRIX.md).
//!
//! Note on determinism: for `standard = 2`, upstream never assigns TPRARR,
//! so the turning-point parity read by GEOKERR is undefined memory. The
//! Rust port sets TPR = 0 for that path and documents it; validation
//! fixtures confirm the standard=2 reference problems are insensitive to
//! this value (see docs/VALIDATION_PLAN.md).

use crate::geokerr::special::gauleg;

/// Arguments for the camera/geodesic initialization.
#[derive(Clone, Debug)]
pub struct CameraArgs {
    pub standard: i32,
    pub a: f64,
    pub mu0: f64,
    pub uout: f64,
    pub u0: f64,
    pub offset: f64,
    pub nro: usize,
    pub nphi: usize,
    pub nup: usize,
}

/// Per-pixel camera data (arrays ALARR..TPRARR of upstream).
#[derive(Clone, Copy, Debug, Default)]
pub struct CameraPixel {
    pub alpha: f64,
    pub beta: f64,
    pub q2: f64,
    pub l: f64,
    pub uf: f64,
    pub muf: f64,
    pub su: f64,
    pub sm: f64,
    pub tpm: i32,
    pub tpr: i32,
}

/// Upstream `initialize_camera_geokerr`. Returns `(u0, offset, pixels)`
/// with pixels in upstream order (`i` = alpha index varies slowest).
#[allow(clippy::too_many_arguments)]
pub fn initialize_camera_geokerr(
    standard: i32,
    a1: f64,
    a2: f64,
    b1: f64,
    b2: f64,
    rcut: f64,
    nrotype: i32,
    nro: usize,
    nphi: usize,
    nup: usize,
    uout: f64,
    mu0: f64,
    a: f64,
) -> (f64, f64, Vec<CameraPixel>) {
    const FAC: f64 = 100.0;
    const INPUT_LIMIT: f64 = 1e-5;
    let zero = 0.0f64;
    let one = 1.0f64;
    let two = 2.0f64;
    let pi = std::f64::consts::PI;

    assert!(
        standard == 1 || standard == 2,
        "initialize_camera_geokerr: only standard=1,2 are used by GRTRANS (got {standard})"
    );
    let _ = uout; // upstream argument, unused in the standard paths

    let ngeo = nro * nphi;

    // effective spin: upstream sets A=0 when |A| < INPUT_LIMIT
    let a_eff = if a.abs() < INPUT_LIMIT { zero } else { a };

    let abmax: f64;
    let ro: Vec<f64>;
    if nrotype == 1 {
        // circular grid: radii are geometric in index
        abmax = rcut.powi(2);
        let r1 = a1;
        let r2 = two.ln();
        // upstream computes Gauss-Legendre nodes then overwrites them
        let (_nodes, _weights) = gauleg(r1, r2, nro);
        ro = (1..=nro)
            .map(|i| {
                // exponent computed in single precision as upstream
                // FLOAT(I)/NRO
                let expo = (i as f32) / (nro as f32);
                r1 * (rcut / r1).powf(expo as f64)
            })
            .collect();
    } else {
        // rectangular grid
        abmax = abmax_rect(a1, a2, b1, b2);
        ro = Vec::new();
    }

    let uplus = one / (one + (one - a_eff * a_eff).sqrt());
    let _umini = one - (one - a_eff * a_eff).sqrt();
    // keep 1/U0 finite
    let u0 = (1e-4f64).min(one / (FAC * abmax));

    let mut offset = 0.5f64;
    if nup == 1 {
        offset = 1e-8;
    }

    let mut pixels = Vec::with_capacity(ngeo);
    if standard == 2 {
        for ii in 1..=ngeo {
            let i = (ii - 1) / nphi + 1;
            let j = (ii - 1) % nphi + 1;
            let (alpha, beta);
            if nrotype > 1 {
                alpha = a1 + (a2 - a1) * (i as f64 - 1.0 + 0.5) / nro as f64;
                beta = b1 + (b2 - b1) * (j as f64 - 1.0 + 0.5) / nphi as f64;
            } else {
                let phi = if nphi != 1 {
                    two * pi * (j as f64 - 1.0 + 0.5) / nphi as f64
                } else {
                    0.0
                };
                alpha = ro[i - 1] * phi.cos();
                beta = ro[i - 1] * phi.sin();
            }
            let l = -alpha * (one - mu0 * mu0).sqrt();
            let mut q2 = beta * beta - (a_eff * a_eff - alpha * alpha) * mu0 * mu0;
            if q2.abs() < INPUT_LIMIT.powi(2) {
                q2 = zero;
            }
            let mut l = l;
            if l.abs() < INPUT_LIMIT {
                l = zero;
            }
            let su = one;
            let sm = if beta > zero && mu0 < one { one } else { -one };
            let tpm = ((if mu0 >= 0.0 { 1.0 } else { -1.0 }) * sm + one) / two;
            let muf = 0.0f64;
            if nup == 1 {
                offset = 0.0;
            }
            pixels.push(CameraPixel {
                alpha,
                beta,
                q2,
                l,
                uf: uplus,
                muf,
                su,
                sm,
                tpm: tpm as i32,
                // upstream leaves TPRARR unassigned here; we use 0
                tpr: 0,
            });
        }
    } else {
        // standard = 1: solve for muf given u0, uf and mu0
        for ii in 1..=ngeo {
            let i = (ii - 1) / nphi + 1;
            let j = (ii - 1) % nphi + 1;
            let uf = uplus;
            let tpr = 1i32;
            let (alpha, beta);
            if nrotype > 1 {
                alpha = a1 + (a2 - a1) * (i as f64 - 1.0 + 0.5) / nro as f64;
                beta = b1 + (b2 - b1) * (j as f64 - 1.0 + 0.5) / nphi as f64;
            } else {
                let phi = if nphi != 1 {
                    two * pi * (j as f64 - 1.0 + 0.5) / nphi as f64
                } else {
                    0.0
                };
                alpha = ro[i - 1] * phi.cos();
                beta = ro[i - 1] * phi.sin();
            }
            let l = -alpha * (one - mu0 * mu0).sqrt();
            let mut q2 = beta * beta - (a_eff * a_eff - alpha * alpha) * mu0 * mu0;
            let su = one;
            let sm = if mu0 < one && beta >= zero { one } else { -one };
            if q2.abs() < INPUT_LIMIT.powi(2) {
                q2 = zero;
            }
            let mut l = l;
            if l.abs() < INPUT_LIMIT {
                l = zero;
            }
            let tpm = 0i32;
            pixels.push(CameraPixel {
                alpha,
                beta,
                q2,
                l,
                uf,
                // upstream does not assign MUFARR here; value unused
                muf: 0.0,
                su,
                sm,
                tpm,
                tpr,
            });
        }
    }

    (u0, offset, pixels)
}

/// Upstream `ABMAX=MAX(A1*A1,A2*A2)**2+MAX(B1*B1,B2*B2)**2` for the
/// rectangular grid (note the fourth power on the alpha term).
#[inline]
fn abmax_rect(a1: f64, a2: f64, b1: f64, b2: f64) -> f64 {
    let ma = (a1 * a1).max(a2 * a2);
    let mb = (b1 * b1).max(b2 * b2);
    ma * ma + mb * mb
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thindisk_camera_geometry() {
        // reproduce the THINDISK test camera: standard=2, 100x100, spin 0.9
        let (u0, offset, pixels) = initialize_camera_geokerr(
            2, -21.0, 21.0, -21.0, 21.0, 1.0, 2, 100, 100, 1, 0.01, 0.26, 0.9,
        );
        assert_eq!(pixels.len(), 10000);
        assert_eq!(offset, 0.0);
        // first pixel: alpha=-21+(0.5)/100*42 = -20.79, beta likewise
        assert!((pixels[0].alpha + 20.79).abs() < 1e-12);
        assert!((pixels[0].beta + 20.79).abs() < 1e-12);
        // u0 = min(1e-4, 1/(100*abmax)) with abmax=2*21^4
        let abmax = 21.0f64.powi(4) * 2.0;
        assert!((u0 - (1.0 / (100.0 * abmax))).abs() < 1e-30);
        // l and q2 for the first pixel
        let l_expected = -pixels[0].alpha * (1.0 - 0.26f64 * 0.26).sqrt();
        assert!((pixels[0].l - l_expected).abs() < 1e-15);
    }
}
