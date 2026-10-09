//! Mu-motion integrals I1MU/I2MU/I3MU (upstream `calcimuasym`,
//! `calcimuasymf`, `calcimusym`, `calcimusymf`).
//!
//! `ImuVals` maps the upstream output arguments: `imu0` is IMU0/IMUF,
//! `imum` is IMUM, `rf0` is RF0/RFF and `rfc` is RFC. The `*f` variants take
//! the incoming IMUM and pass it through unchanged (upstream never modifies
//! it in those routines).

use crate::geokerr::elliptic::rf;

/// Outputs of the `calcimu*` routines.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImuVals {
    pub imu0: f64,
    pub imum: f64,
    pub rf0: f64,
    pub rfc: f64,
}

/// Upstream `calcimuasym`: asymmetric roots, starting from MU0.
pub fn calcimuasym(a: f64, mneg: f64, mpos: f64, mu0: f64, muplus: f64) -> ImuVals {
    let fac = a.abs() * muplus;
    let m1 = mneg / mpos;
    let phi0 = ((mpos - mu0 * mu0) / (mpos - mneg)).sqrt().asin();
    let s0 = ((mpos - mu0 * mu0) / (mpos - mneg)).sqrt();
    let ak = (1.0 - m1).sqrt();
    let q0 = (1.0 - s0 * ak) * (1.0 + s0 * ak);
    let rf0 = rf(1.0 - s0 * s0, q0, 1.0);
    let rfc = rf(0.0, m1, 1.0);
    let mut imum = 2.0 * rfc;
    let mut imu0 = s0 * rf0;
    if phi0.tan() < 0.0 {
        imu0 = imum - imu0;
    }
    imu0 /= fac;
    imum = imum / fac / 2.0;
    ImuVals {
        imu0,
        imum,
        rf0,
        rfc,
    }
}

/// Upstream `calcimuasymf`: asymmetric roots, ending at MUF.
pub fn calcimuasymf(a: f64, mneg: f64, mpos: f64, muf: f64, muplus: f64, imum_in: f64) -> ImuVals {
    let fac = a.abs() * muplus;
    let m1 = mneg / mpos;
    let sf = ((mpos - muf * muf) / (mpos - mneg)).sqrt();
    let ak = (1.0 - m1).sqrt();
    let qf = (1.0 - sf * ak) * (1.0 + sf * ak);
    let rff = rf(1.0 - sf * sf, qf, 1.0);
    let mut imuf = sf * rff;
    imuf /= fac;
    ImuVals {
        imu0: imuf,
        imum: imum_in,
        rf0: rff,
        rfc: 0.0,
    }
}

/// Upstream `calcimusym`: symmetric roots, starting from MU0.
pub fn calcimusym(a: f64, mneg: f64, mpos: f64, mu0: f64, muplus: f64) -> ImuVals {
    let fac = a.abs() * (mpos - mneg).sqrt();
    let m1 = -mneg / (mpos - mneg);
    let phi0 = (mu0 / muplus).acos();
    let s0 = phi0.sin();
    let ak = (1.0 - m1).sqrt();
    let q0 = (1.0 - s0 * ak) * (1.0 + s0 * ak);
    let rf0 = rf(1.0 - s0 * s0, q0, 1.0);
    let rfc = rf(0.0, m1, 1.0);
    let mut imum = 2.0 * rfc;
    let mut imu0 = s0 * rf0;
    if phi0.tan() < 0.0 {
        imu0 = imum - imu0;
    }
    imu0 /= fac;
    imum /= fac;
    ImuVals {
        imu0,
        imum,
        rf0,
        rfc,
    }
}

/// Upstream `calcimusymf`: symmetric roots, ending at MUF.
pub fn calcimusymf(a: f64, mneg: f64, mpos: f64, muf: f64, muplus: f64, imum_in: f64) -> ImuVals {
    let fac = a.abs() * (mpos - mneg).sqrt();
    let m1 = -mneg / (mpos - mneg);
    let phif = (muf / muplus).acos();
    let sf = phif.sin();
    let ak = (1.0 - m1).sqrt();
    let qf = (1.0 - sf * ak) * (1.0 + sf * ak);
    let rff = rf(1.0 - sf * sf, qf, 1.0);
    let mut imuf = sf * rff;
    if phif.tan() < 0.0 {
        imuf = imum_in * fac - imuf;
    }
    imuf /= fac;
    ImuVals {
        imu0: imuf,
        imum: imum_in,
        rf0: rff,
        rfc: 0.0,
    }
}
