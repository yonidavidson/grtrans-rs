//! Relativistic Maxwell-Jüttner (`calc_maxjutt.f90`) and Maxwell-Compton
//! (`calc_maxcomp.f90`) electron-distribution decompositions.
//!
//! Both are weighted sums of the thermal polarized synchrotron
//! (`polsynchth`, already validated) over a geometric temperature ladder.

use crate::polsynch::polsynchth;

/// Shared ladder construction: normalized weights, `delta^i`, `dwsum` and the
/// per-point `tcgs_min = tcgs / dwsum`.
fn ladder(weights_raw: &[f64], delta: f64, tcgs: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let sum: f64 = weights_raw.iter().sum();
    let weights: Vec<f64> = weights_raw.iter().map(|w| w / sum).collect();
    let n = weights.len();
    let delta_arr: Vec<f64> = (0..n).map(|i| delta.powf(i as f64)).collect();
    let dwsum: f64 = weights
        .iter()
        .zip(delta_arr.iter())
        .map(|(w, d)| w * d)
        .sum();
    let tcgs_min: Vec<f64> = tcgs.iter().map(|t| t / dwsum).collect();
    (weights, delta_arr, tcgs_min)
}

fn weighted_sum(
    nu: &[f64],
    ncgs: &[f64],
    bcgs: &[f64],
    incang: &[f64],
    weights: &[f64],
    delta_arr: &[f64],
    tcgs_min: &[f64],
) -> Vec<[f64; 11]> {
    let np = nu.len();
    let mut ktotal = vec![[0.0f64; 11]; np];
    for (w, d) in weights.iter().zip(delta_arr.iter()) {
        let n_i: Vec<f64> = ncgs.iter().map(|v| w * v).collect();
        let t_i: Vec<f64> = tcgs_min.iter().map(|v| d * v).collect();
        let ktemp = polsynchth(nu, &n_i, bcgs, &t_i, incang);
        for (j, row) in ktemp.iter().enumerate() {
            for c in 0..11 {
                ktotal[j][c] += row[c];
            }
        }
    }
    ktotal
}

/// Upstream `calc_maxjutt_subroutine`. `args = [delta, w1, ..., wN]`.
pub fn calc_maxjutt(
    nu: &[f64],
    ncgs: &[f64],
    bcgs: &[f64],
    tcgs: &[f64],
    incang: &[f64],
    args: &[f64],
) -> Vec<[f64; 11]> {
    let delta = args[0];
    let (weights, delta_arr, tcgs_min) = ladder(&args[1..], delta, tcgs);
    weighted_sum(nu, ncgs, bcgs, incang, &weights, &delta_arr, &tcgs_min)
}

/// Upstream `calc_maxcomp_subroutine`. `args = [delta, selection, w1, ..., wN]`.
/// If `selection` (as a Fortran default-real, i.e. f32, truncated to integer)
/// is in `1..=N`, the emission columns are replaced by that single
/// Maxwellian while the total absorption is kept.
pub fn calc_maxcomp(
    nu: &[f64],
    ncgs: &[f64],
    bcgs: &[f64],
    tcgs: &[f64],
    incang: &[f64],
    args: &[f64],
) -> Vec<[f64; 11]> {
    let delta = args[0];
    // upstream `selection` is a default-real (f32) variable; INT truncates.
    let isel = (args[1] as f32) as i32;
    let (weights, delta_arr, tcgs_min) = ladder(&args[2..], delta, tcgs);
    let mut ktotal = weighted_sum(nu, ncgs, bcgs, incang, &weights, &delta_arr, &tcgs_min);
    let n = weights.len();
    if isel > 0 && (isel as usize) <= n {
        let idx = (isel - 1) as usize;
        let (w, d) = (weights[idx], delta_arr[idx]);
        let n_i: Vec<f64> = ncgs.iter().map(|v| w * v).collect();
        let t_i: Vec<f64> = tcgs_min.iter().map(|v| d * v).collect();
        let ktemp = polsynchth(nu, &n_i, bcgs, &t_i, incang);
        for (j, row) in ktemp.iter().enumerate() {
            for c in 0..4 {
                ktotal[j][c] = row[c];
            }
        }
    }
    ktotal
}
