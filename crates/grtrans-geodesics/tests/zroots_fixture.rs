//! ZROOTS validation against the Fortran fixture.
use grtrans_geodesics::geokerr::special::zroots;
use num_complex::Complex64;

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference/fixtures/fortran/test_zroots.txt"
);

#[test]
fn zroots_quartic_matches_fortran() {
    let text = std::fs::read_to_string(FIXTURE).expect("test_zroots fixture missing");
    let mut cc = 0.0;
    let mut dd = 0.0;
    let mut ee = 0.0;
    let mut quartic = Vec::new();
    let mut in_quartic = false;
    for line in text.lines() {
        if line.starts_with("# cc dd ee") {
            let vals: Vec<f64> = line
                .split_whitespace()
                .skip(4)
                .map(|t| t.parse().unwrap())
                .collect();
            cc = vals[0];
            dd = vals[1];
            ee = vals[2];
        } else if line.starts_with("# quartic") {
            in_quartic = true;
        } else if line.starts_with('#') {
            in_quartic = false;
        } else if in_quartic {
            let vals: Vec<f64> = line
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            quartic.push((vals[0], vals[1]));
        }
    }
    assert_eq!(quartic.len(), 4);

    let c: [Complex64; 5] = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(cc, 0.0),
        Complex64::new(dd, 0.0),
        Complex64::new(ee, 0.0),
    ];
    let roots = zroots(&c, 4, true);
    assert_eq!(roots.len(), 4);
    for (i, (re, im)) in quartic.iter().enumerate() {
        let r = roots[i];
        let scale = re.abs().max(1e-300);
        assert!(
            (r.re - re).abs() / scale < 1e-13,
            "root {i} real: rust {} vs fortran {re}",
            r.re
        );
        assert!(
            (r.im - im).abs() < 1e-13,
            "root {i} imag: rust {} vs {im}",
            r.im
        );
    }
}
