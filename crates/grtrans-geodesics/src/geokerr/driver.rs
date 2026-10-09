//! The GEOKERR driver: computes NUP points along one null geodesic.
//!
//! Direct translation of subroutine GEOKERR in `geokerr_wrapper.f`
//! (lines 357-604). Upstream passes scalars by reference and modifies
//! several of them (NUP, UF, TPR, TPM, MUF); the Rust signature returns the
//! updated values in [`GeokerrResult`].

use crate::geokerr::elliptic::asech;
use crate::geokerr::geomu::{findmroots, geomu, indep_muf, GeomuInputs};
use crate::geokerr::phitime::{geophitime, PhitimeState};
use crate::geokerr::radial::{geor, GeorInputs};

/// Result of one geokerr call. `ufi`, `mufi`, `dti`, `dphi`, `tpmi`, `tpri`
/// and `lambdai` have length `nup_alloc = nup_in + npts_extra - 2*kext`,
/// matching the upstream output arrays; only the first `nup` entries are
/// meaningful (upstream's modified NUP).
#[derive(Clone, Debug, Default)]
pub struct GeokerrResult {
    /// modified NUP (upstream updates the argument in place)
    pub nup: usize,
    /// possibly modified UF, TPR, TPM, MUF
    pub uf: f64,
    pub tpr: i32,
    pub tpm: i32,
    pub muf: f64,
    pub ufi: Vec<f64>,
    pub mufi: Vec<f64>,
    pub dti: Vec<f64>,
    pub dphi: Vec<f64>,
    pub tpmi: Vec<i32>,
    pub tpri: Vec<i32>,
    pub lambdai: Vec<f64>,
}

/// Fortran `SIGN(1d0, x)`.
#[inline]
fn fsign(x: f64) -> f64 {
    if x >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// Fortran `INT(x)`.
#[inline]
fn f_int(x: f64) -> i64 {
    x as i64
}

/// Upstream `geokerr`. `kext`/`npts_extra` mirror upstream's KEXT/NPTS
/// (the MUFILL extension parameters from `initialize_pixels`).
#[allow(clippy::too_many_arguments)]
pub fn geokerr(
    u0: f64,
    uf_in: f64,
    uout: f64,
    mu0: f64,
    muf_in: f64,
    a: f64,
    l: f64,
    q2: f64,
    alpha: f64,
    beta: f64,
    tpm_in: i32,
    tpr_in: i32,
    su: f64,
    sm: f64,
    nup_in: usize,
    offset: f64,
    phit: bool,
    usegeor: bool,
    mufill: bool,
    kext: usize,
    npts_extra: usize,
) -> GeokerrResult {
    let one = 1.0f64;
    let two = 2.0f64;
    let l2 = l * l;
    let uplus = one / (one + (one - a * a).sqrt());

    let nup_alloc = nup_in + npts_extra - 2 * kext;
    let mut out = GeokerrResult {
        ufi: vec![0.0; nup_alloc],
        mufi: vec![0.0; nup_alloc],
        dti: vec![0.0; nup_alloc],
        dphi: vec![0.0; nup_alloc],
        tpmi: vec![0; nup_alloc],
        tpri: vec![0; nup_alloc],
        lambdai: vec![0.0; nup_alloc],
        ..Default::default()
    };

    // shared mutable state (upstream's scalars)
    let mut uf = uf_in;
    let mut muf = muf_in;
    let mut tpm = tpm_in;
    let mut tpr = tpr_in;
    let mut npts = npts_extra;
    let mut firstpt = false;

    // cached first-point data shared between GEOMU / GEOR / GEOPHITIME
    let mut gmu_cache = GeomuInputs::first();
    let mut gr_cache = GeorInputs::default();
    let mut phi_state = PhitimeState::default();
    // rffmu2 is computed fresh by GEOR/GEOMU whenever PHIT is true; it is
    // not part of the cached first-point state (matching upstream's
    // argument usage).
    let mut rffmu2 = 0.0f64;

    if usegeor {
        // Solve for uf at equal steps between mu0 and muf by calling GEOR
        // assuming no mu turning points are present.
        let mut k = 1usize;
        let mun = mu0 + (k as f64 - offset) * (muf - mu0) / nup_in as f64;
        let mut lambda = 0.0f64;
        let mut tmu = 0.0f64;
        let mut tu = 0.0f64;
        let mut phimu = 0.0f64;
        let mut phiu = 0.0f64;
        if mu0 != 0.0 || beta != 0.0 {
            let res = geor(
                u0, uf, mu0, mun, a, l, l2, q2, tpm, tpr, su, sm, &gr_cache, phit, true,
            );
            uf = res.uf;
            tpr = res.tpr;
            rffmu2 = res.rffmu2;
            gr_cache = res.cache();
            if phit {
                let p = geophitime(
                    u0, uf, mu0, mun, a, l, l2, q2, tpm, tpr, su, sm, res.iu, res.h1, res.u1,
                    res.u2, res.u3, res.u4, res.rffu0, res.rffu1, res.rffmu1, rffmu2, res.rffmu3,
                    res.ncase, &phi_state, true,
                );
                phimu = p.phimu;
                tmu = p.tmu;
                phiu = p.phiu;
                tu = p.tu;
                lambda = p.lambda;
                phi_state = p.state;
            }
        } else {
            uf = u0;
        }
        out.ufi[0] = uf;
        out.mufi[0] = mun;
        out.dti[0] = tmu + tu;
        out.dphi[0] = phimu + phiu;
        out.tpmi[0] = tpm;
        out.tpri[0] = tpr;
        out.lambdai[0] = lambda;

        for k in 2..=nup_in {
            let mun = mu0 + (k as f64 - offset) * (muf - mu0) / nup_in as f64;
            tpm = ((fsign(mu0) * sm + one) / two) as i32;
            let (dti, dphi, lambda_k, tpri_k);
            if mu0 != 0.0 || beta != 0.0 {
                let res = geor(
                    u0, uf, mu0, mun, a, l, l2, q2, tpm, tpr, su, sm, &gr_cache, phit, false,
                );
                uf = res.uf;
                tpr = res.tpr;
                rffmu2 = res.rffmu2;
                gr_cache.rffu0 = res.rffu0;
                gr_cache.iu0 = res.iu0;
                if phit {
                    let p = geophitime(
                        u0, uf, mu0, mun, a, l, l2, q2, tpm, tpr, su, sm, res.iu, res.h1, res.u1,
                        res.u2, res.u3, res.u4, res.rffu0, res.rffu1, res.rffmu1, rffmu2,
                        res.rffmu3, res.ncase, &phi_state, true,
                    );
                    phimu = p.phimu;
                    tmu = p.tmu;
                    phiu = p.phiu;
                    tu = p.tu;
                    lambda = p.lambda;
                    phi_state = p.state;
                }
                dti = tmu + tu;
                dphi = phimu + phiu;
                lambda_k = lambda;
                tpri_k = tpr;
            } else {
                uf = u0;
                dti = 0.0;
                dphi = 0.0;
                lambda_k = 0.0;
                tpri_k = tpr;
            }
            out.ufi[k - 1] = uf;
            out.mufi[k - 1] = mun;
            out.dti[k - 1] = dti;
            out.dphi[k - 1] = dphi;
            out.tpmi[k - 1] = tpm;
            out.tpri[k - 1] = tpri_k;
            out.lambdai[k - 1] = lambda_k;
        }
        out.nup = nup_in;
        out.uf = uf;
        out.muf = muf;
        out.tpr = tpr;
        out.tpm = tpm;
        return out;
    }

    // ------------------------------------------------------------------
    // USEGEOR = false: integrate from u0 outwards
    // ------------------------------------------------------------------
    let mut nup = nup_in;
    let first = geomu(
        u0, uf, mu0, muf, a, l, l2, q2, tpm, tpr, su, sm, &gmu_cache, phit, true,
    );
    uf = first.uf;
    tpr = first.tpr;
    tpm = first.tpm;
    muf = first.muf;
    gmu_cache = first.cache;
    gr_cache = GeorInputs {
        h1: first.h1,
        u1: first.u1,
        u2: first.u2,
        u3: first.u3,
        u4: first.u4,
        rffu0: first.rffu0,
        rffu1: first.rffu1,
        rffmu1: first.rffmu1,
        rffmu2: first.rffmu2,
        rffmu3: first.rffmu3,
        iu0: first.iu0,
        i1mu: first.i1mu,
        i2mu: 0.0,
        i3mu: first.i3mu,
        ncase: first.ncase,
    };
    let ncase = first.ncase;

    if (ncase > 2 && ncase < 7)
        || (ncase == 1 && su < 0.0)
        || (ncase == 7 && su < 0.0)
        || (ncase == 2 && su > 0.0)
        || (ncase == 8 && su > 0.0)
    {
        tpr = 0;
    }

    let ub;
    if ncase > 2 && ncase < 7 {
        ub = if su > 0.0 { uplus } else { 0.0 };
    } else if ncase == 1 || ncase == 7 {
        ub = first.u2;
        if tpr == 1 && uf == ub {
            uf = uout;
        }
    } else {
        ub = first.u3;
        if tpr == 1 && uf == ub {
            uf = uout;
        }
    }

    let tpr1 = 0i32;
    let mut du = fsign(ub - uout) * su * ((ub - uout) + (two * tpr as f64 - one) * (ub - uf));
    du /= nup as f64;
    // Fortran: KMAX = MIN(INT((UB-UOUT)/DU+OFFSET), NUP); INT truncates.
    // The clamp at zero only affects pathological inputs where upstream
    // would index out of bounds.
    let mut kmax = (f_int((ub - uout) / du + offset)).min(nup as i64).max(0);
    if du == 0.0 {
        kmax = 0;
    }
    let kmax = kmax as usize;

    if fsign(uout - ub) != fsign(uout - u0) || uout == u0 {
        if kmax != 0 {
            // First point computed separately so that the once-per-geodesic
            // integrals are computed only once.
            let k = 1usize;
            let un = uout + (k as f64 - offset) * du;
            let mun = 0.0f64; // see note below: upstream passes MUN here, which
                              // GEOMU overwrites with the computed mu
            let res = geomu(
                u0, un, mu0, mun, a, l, l2, q2, tpm, tpr1, su, sm, &gmu_cache, phit, false,
            );
            gmu_cache = res.cache;
            rffmu2 = res.rffmu2;
            let mut lambda = 0.0;
            let mut dti = 0.0;
            let mut dphi = 0.0;
            if phit {
                let p = geophitime(
                    u0, un, mu0, res.muf, a, l, l2, q2, res.tpm, tpr1, su, sm, res.iu, res.h1,
                    res.u1, res.u2, res.u3, res.u4, res.rffu0, res.rffu1, res.rffmu1, rffmu2,
                    res.rffmu3, res.ncase, &phi_state, true,
                );
                dti = p.tmu + p.tu;
                dphi = p.phimu + p.phiu;
                lambda = p.lambda;
                phi_state = p.state;
            }
            out.ufi[0] = un;
            out.mufi[0] = res.muf;
            out.dti[0] = dti;
            out.dphi[0] = dphi;
            out.tpmi[0] = res.tpm;
            out.tpri[0] = 0;
            out.lambdai[0] = lambda;
        }
        // Trace from u0 to uf or turning point
        for k in 2..=kmax {
            let un = uout + (k as f64 - offset) * du;
            let res = geomu(
                u0, un, mu0, 0.0, a, l, l2, q2, tpm, tpr1, su, sm, &gmu_cache, phit, false,
            );
            gmu_cache = res.cache;
            rffmu2 = res.rffmu2;
            let mut lambda = 0.0;
            let mut dti = 0.0;
            let mut dphi = 0.0;
            if phit {
                let p = geophitime(
                    u0, un, mu0, res.muf, a, l, l2, q2, res.tpm, tpr1, su, sm, res.iu, res.h1,
                    res.u1, res.u2, res.u3, res.u4, res.rffu0, res.rffu1, res.rffmu1, rffmu2,
                    res.rffmu3, res.ncase, &phi_state, false,
                );
                dti = p.tmu + p.tu;
                dphi = p.phimu + p.phiu;
                lambda = p.lambda;
            }
            out.ufi[k - 1] = un;
            out.mufi[k - 1] = res.muf;
            out.dti[k - 1] = dti;
            out.dphi[k - 1] = dphi;
            out.tpmi[k - 1] = res.tpm;
            out.tpri[k - 1] = 0;
            out.lambdai[k - 1] = lambda;
        }
        if tpr == 1 && kmax > kext {
            if nup == 1 {
                firstpt = true;
            }
            if mufill {
                // Fill points near the u turning point using mu as the
                // independent variable.
                let k = kmax + kext + 1;
                let un = two * ub - (uout + (k as f64 - offset) * du);
                let res = geomu(
                    uf, un, mu0, 0.0, a, l, l2, q2, tpm, tpr, su, sm, &gmu_cache, phit, false,
                );
                gmu_cache = res.cache;
                tpr = res.tpr;
                let tpm1 = res.tpm;
                muf = res.muf;
                let (muminus, muplus) = findmroots(q2, l, a, mu0);
                for kk in 1..=npts {
                    let (mun, tpmk) = indep_muf(
                        muminus,
                        muplus,
                        out.mufi[kmax - kext - 1],
                        muf,
                        out.tpmi[kmax - kext - 1],
                        tpm1,
                        sm,
                        kk as i32,
                        npts as i32,
                        offset,
                    );
                    tpm = tpmk;
                    let res2 = geor(
                        u0, un, mu0, mun, a, l, l2, q2, tpm, tpr, su, sm, &gr_cache, phit, false,
                    );
                    tpr = res2.tpr;
                    rffmu2 = res2.rffmu2;
                    // GEOR's UF argument is in/out: upstream updates UN in
                    // place, and the next fill iteration uses the new value.
                    let un = res2.uf;
                    let mut lambda = 0.0;
                    let mut dti = 0.0;
                    let mut dphi = 0.0;
                    if phit {
                        let p = geophitime(
                            u0,
                            un,
                            mu0,
                            mun,
                            a,
                            l,
                            l2,
                            q2,
                            tpm,
                            tpr,
                            su,
                            sm,
                            res2.iu,
                            res2.h1,
                            res2.u1,
                            res2.u2,
                            res2.u3,
                            res2.u4,
                            res2.rffu0,
                            res2.rffu1,
                            res2.rffmu1,
                            rffmu2,
                            res2.rffmu3,
                            res2.ncase,
                            &phi_state,
                            false,
                        );
                        dti = p.tmu + p.tu;
                        dphi = p.phimu + p.phiu;
                        lambda = p.lambda;
                    }
                    let idx = kk + kmax - kext - 1;
                    // Upstream writes up to KEXT elements past the official
                    // array length here (the values are never read because
                    // the returned NUP excludes them). Guard instead of
                    // reproducing the out-of-bounds heap write.
                    if idx < out.ufi.len() {
                        out.ufi[idx] = un;
                        out.mufi[idx] = mun;
                        out.dti[idx] = dti;
                        out.dphi[idx] = dphi;
                        out.tpmi[idx] = tpm;
                        out.tpri[idx] = tpr;
                        out.lambdai[idx] = lambda;
                    }
                }
                npts -= 2 * kext;
            }
            // Trace from the turning point back to uf
            for k in (kmax + 1 + kext)..=nup {
                let un = two * ub - (uout + (k as f64 - offset) * du);
                let res = geomu(
                    uf, un, mu0, 0.0, a, l, l2, q2, tpm, tpr, su, sm, &gmu_cache, phit, false,
                );
                gmu_cache = res.cache;
                tpr = res.tpr;
                rffmu2 = res.rffmu2;
                let mut lambda = 0.0;
                let mut dti = 0.0;
                let mut dphi = 0.0;
                if phit {
                    let p = geophitime(
                        uf, un, mu0, res.muf, a, l, l2, q2, res.tpm, tpr, su, sm, res.iu, res.h1,
                        res.u1, res.u2, res.u3, res.u4, res.rffu0, res.rffu1, res.rffmu1, rffmu2,
                        res.rffmu3, res.ncase, &phi_state, firstpt,
                    );
                    dti = p.tmu + p.tu;
                    dphi = p.phimu + p.phiu;
                    lambda = p.lambda;
                }
                let idx = k + npts - 1;
                out.ufi[idx] = un;
                out.mufi[idx] = res.muf;
                out.dti[idx] = dti;
                out.dphi[idx] = dphi;
                out.tpmi[idx] = res.tpm;
                out.tpri[idx] = 1;
                out.lambdai[idx] = lambda;
            }
            nup += npts;
        }
    } else {
        nup = 0;
    }

    out.nup = nup;
    out.uf = uf;
    out.muf = muf;
    out.tpr = tpr;
    out.tpm = tpm;
    out
}

/// Keep the `asech` import used by future extensions (GEOR uses it).
#[allow(dead_code)]
fn _unused(x: f64) -> f64 {
    asech(x)
}
