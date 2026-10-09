//! Chandrasekhar (1960) table XXIV: degree of polarization of radiation
//! emerging from a semi-infinite electron-scattering atmosphere.
//!
//! Direct translation of `chandra_tab24.f90` (upstream GRTRANS). The table
//! data (`ch24_vals.txt`) is embedded from the upstream repository and is
//! loaded once on first use.

use crate::interpolate::get_weight;
use std::sync::OnceLock;

/// Upstream table data, verbatim from `ch24_vals.txt`.
pub const CH24_VALS: &str = include_str!("../data/ch24_vals.txt");

#[derive(Debug)]
pub struct ChandraTable {
    /// mu = cos(viewing angle), ascending
    pub mu: Vec<f64>,
    /// total intensity
    pub i: Vec<f64>,
    /// degree of polarization
    pub delta: Vec<f64>,
}

impl ChandraTable {
    /// Parse the upstream table format: `npts` on the first line (the rest
    /// of that line is discarded, as with Fortran list-directed input),
    /// followed by three arrays: `mu`, `I`, `delta`.
    pub fn parse(data: &str) -> Self {
        let lines: Vec<&str> = data.lines().filter(|l| !l.trim().is_empty()).collect();
        let first: Vec<&str> = lines[0].split_whitespace().collect();
        let npts: usize = first[0].parse().expect("ch24 npts");
        let mut values: Vec<f64> = Vec::with_capacity(npts * 3);
        for line in lines.iter().skip(1) {
            for tok in line.split_whitespace() {
                values.push(tok.parse().expect("ch24 value"));
            }
        }
        assert!(values.len() >= npts * 3, "ch24 table truncated");
        let mu = values[..npts].to_vec();
        let i = values[npts..2 * npts].to_vec();
        let delta = values[2 * npts..3 * npts].to_vec();
        ChandraTable { mu, i, delta }
    }
}

/// The embedded table (loaded once).
pub fn table() -> &'static ChandraTable {
    static TABLE: OnceLock<ChandraTable> = OnceLock::new();
    TABLE.get_or_init(|| ChandraTable::parse(CH24_VALS))
}

/// Interpolate intensity and polarization degree at `mu`
/// (upstream `interp_chandra_tab24`). Returned as `(I, delta)`.
pub fn interp_chandra_tab24(mu: f64) -> (f64, f64) {
    let t = table();
    let (weight, j) = get_weight(&t.mu, mu, 0);
    let j0 = (j - 1) as usize;
    let i = (1.0 - weight) * t.i[j0] + weight * t.i[j0 + 1];
    let del = (1.0 - weight) * t.delta[j0] + weight * t.delta[j0 + 1];
    (i, del)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_parses_and_interpolates() {
        let t = table();
        assert_eq!(t.mu.len(), 21);
        assert_eq!(t.mu[0], 0.0);
        assert_eq!(t.mu[20], 1.0);
        // endpoints exact
        let (i0, d0) = interp_chandra_tab24(0.0);
        assert!((i0 - 0.414410).abs() < 1e-6);
        assert!((d0 - 0.117130).abs() < 1e-6);
        let (i1, d1) = interp_chandra_tab24(1.0);
        assert!((i1 - 1.26938).abs() < 1e-5);
        assert!((d1 - 0.0).abs() < 1e-6);
        // midpoint
        let (i5, _) = interp_chandra_tab24(0.5);
        assert!((i5 - 0.86637).abs() < 1e-6);
    }
}

/// The embedded table as single-precision values (upstream reads
/// `ch_mu, ch_I, ch_delta` as default `real`).
fn table_f32() -> &'static ChandraTableF32 {
    static TABLE: OnceLock<ChandraTableF32> = OnceLock::new();
    TABLE.get_or_init(|| {
        let t = table();
        ChandraTableF32 {
            mu: t.mu.iter().map(|v| *v as f32).collect(),
            i: t.i.iter().map(|v| *v as f32).collect(),
            delta: t.delta.iter().map(|v| *v as f32).collect(),
        }
    })
}

/// Single-precision table (upstream `ch_mu`, `ch_I`, `ch_delta`).
#[derive(Debug)]
pub struct ChandraTableF32 {
    pub mu: Vec<f32>,
    pub i: Vec<f32>,
    pub delta: Vec<f32>,
}

/// Interpolate intensity and polarization degree in single precision,
/// exactly as upstream `interp_chandra_tab24` does (its arguments are
/// default `real`).
pub fn interp_chandra_tab24_f32(mu: f32) -> (f32, f32) {
    let t = table_f32();
    let (weight, j) = crate::interpolate::get_weight_f32(&t.mu, mu, 0);
    let j0 = (j - 1) as usize;
    let i = (1.0 - weight) * t.i[j0] + weight * t.i[j0 + 1];
    let del = (1.0 - weight) * t.delta[j0] + weight * t.delta[j0 + 1];
    (i, del)
}
