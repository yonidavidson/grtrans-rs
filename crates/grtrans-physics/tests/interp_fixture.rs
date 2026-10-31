//! INTERP (`interpemis_noabs`) validated against the upstream Fortran
//! routine (`test_interp`). Also documents the below-grid edge deviation:
//! upstream reads `freqarr(0)` (out of bounds) and returns NaN; we clamp to
//! the first segment instead (a finite value).

use grtrans_physics::emissivity::interpemis_noabs;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_interp.txt"
);

#[test]
fn interp_matches_fortran() {
    let nf = 4usize;
    let np = 5usize;
    let freqarr = [1.0e10, 1.0e11, 1.0e12, 1.0e13];
    let mut jnu = vec![0.0f64; np * nf];
    for p in 0..np {
        for f in 0..nf {
            jnu[p * nf + f] = ((p + 1) as f64) * 10f64.powi(f as i32);
        }
    }
    let nu = [5.0e10, 1.0e11, 3.0e12, 1.0e9, 1.0e14];
    let out = interpemis_noabs(&nu, &freqarr, &jnu);
    let expected: Vec<(f64, f64)> = std::fs::read_to_string(FIXTURE)
        .unwrap()
        .lines()
        .map(|l| {
            let t: Vec<f64> = l
                .split_whitespace()
                .map(|x| x.replace(['d', 'D'], "e").parse().unwrap())
                .collect();
            (t[0], t[1])
        })
        .collect();
    assert_eq!(expected.len(), np);
    for (i, (nu_i, k_f)) in expected.iter().enumerate() {
        assert_eq!(*nu_i, nu[i]);
        let ours = out[i][0];
        if k_f.is_nan() {
            // below-grid: upstream is NaN (out-of-bounds read); we clamp to a
            // finite value (documented deviation)
            assert!(ours.is_finite(), "below-grid should be finite, got {ours}");
            println!("below-grid nu={nu_i:e}: ours={ours:e} (upstream NaN, clamped)");
        } else {
            let rel = (ours - k_f).abs() / k_f.abs().max(1e-300);
            println!("nu={nu_i:e}: ours={ours:e} fortran={k_f:e} rel={rel:e}");
            assert!(rel < 1e-5, "nu={nu_i:e}: rel {rel:e} exceeds 1e-5");
        }
    }
}
